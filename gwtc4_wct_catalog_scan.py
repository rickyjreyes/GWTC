#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
GWTC-4 WCT Catalog / Posterior Branch Scanner
---------------------------------------------

Purpose
-------
Scan GWTC-4 posterior-sample HDF5 files for WCT-style log-domain and
integer-branch structure.

This script is designed for the LVK / IGWN GWTC-4 posterior release files named like:

    IGWN-GWTC4p0-...-GW231123_135430-combined_PEDataRelease.hdf5
    IGWN-GWTC4p0-...-PESummaryTable.hdf5

It is NOT a strain-waveform scanner. It is the first, safer catalog/posterior layer:
    - read per-event posterior samples,
    - extract masses/spins if available,
    - reduce each event to medians and widths,
    - scan catalog-level ordered sequences for log-periodic residual structure,
    - test branch ratios such as STAR 6:17 and LHC/Koide 10:15:20.

Core WCT idea tested here
-------------------------
For a positive catalog variable x, define:

    ell = ln(x)

For an ordered catalog sequence y(ell), test:

    y(ell) = baseline(ell) + a cos(k ell) + b sin(k ell)

Then compare best k to integer branch grids:

    k_n = 2*pi*n / DeltaEll

and ratio branches:
    (6, 17)          STAR-transposed branch
    (10, 15, 20)    LHC / Koide Class-I branch
    (20/3,15,40/3) folded Q=4/9 branch

Important
---------
This script is exploratory unless you lock a variable/order/branch rule before
running on a held-out subset.

Install
-------
    pip install numpy pandas scipy h5py matplotlib

Optional GPU:
    This script uses NumPy/SciPy by default because catalog scans are small.
    It does not require CuPy.

Run
---
    # Verify downloads against md5sums.txt
    python gwtc4_wct_catalog_scan.py --manifest md5sums.txt --verify-only

    # Inspect HDF5 paths
    python gwtc4_wct_catalog_scan.py --input "*.hdf5" --inspect-only --inspect-limit 3

    # Full scan
    python gwtc4_wct_catalog_scan.py --input "*combined_PEDataRelease.hdf5" --null-n 5000

    # With manifest verification
    python gwtc4_wct_catalog_scan.py --manifest md5sums.txt --input "*combined_PEDataRelease.hdf5" --null-n 5000

Outputs
-------
    outputs_gwtc4_wct/
        gwtc4_event_summary.csv
        gwtc4_scan_results.csv
        gwtc4_branch_match_results.csv
        gwtc4_ratio_results.csv
        gwtc4_run_summary.json
        scans/*.csv
        nulls/*.csv
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
import os
import re
import sys
import warnings
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import h5py
import numpy as np
import pandas as pd

try:
    from scipy.signal import find_peaks
except Exception:
    find_peaks = None


# =============================================================================
# Defaults
# =============================================================================

OUTDIR_DEFAULT = "outputs_gwtc4_wct"

RNG_SEED_DEFAULT = 12345

K_MIN_DEFAULT = 0.5
K_MAX_DEFAULT = 40.0
N_K_DEFAULT = 3000

BASELINE_DEGREE_DEFAULT = 2
NULL_N_DEFAULT = 5000

MIN_EVENTS_DEFAULT = 12

EPS = 1e-12

# LHC active-domain length from your LHCb paper.
DELTA_ELL_LHC = 4.780150335923678

# Branches tested as ratios / integer templates.
BRANCHES = {
    "star_transposed_6_17": [6.0, 17.0],
    "koide_10_15_20": [10.0, 15.0, 20.0],
    "folded_4over9": [20.0 / 3.0, 15.0, 40.0 / 3.0],
    "integer_1_to_30": list(map(float, range(1, 31))),
    "integer_6_to_17": list(map(float, range(6, 18))),
    "integer_10_to_22": list(map(float, range(10, 23))),
}

# Candidate dataset names. HDF5 releases differ by event and waveform.
# The scanner searches by suffix/name match and by common PE summary names.
PARAM_ALIASES = {
    # masses
    "mass_1_source": [
        "mass_1_source", "m1_source", "srcmass1", "source_mass_1",
        "mass1_source", "mass_1_source_solar_mass",
    ],
    "mass_2_source": [
        "mass_2_source", "m2_source", "srcmass2", "source_mass_2",
        "mass2_source", "mass_2_source_solar_mass",
    ],
    "mass_1_detector": [
        "mass_1", "mass_1_detector", "m1", "mass1", "detector_mass_1",
    ],
    "mass_2_detector": [
        "mass_2", "mass_2_detector", "m2", "mass2", "detector_mass_2",
    ],
    "chirp_mass_source": [
        "chirp_mass_source", "mc_source", "source_chirp_mass", "chirp_mass_source_solar_mass",
    ],
    "chirp_mass_detector": [
        "chirp_mass", "chirp_mass_detector", "mc", "detector_chirp_mass",
    ],
    "total_mass_source": [
        "total_mass_source", "mtotal_source", "source_total_mass",
    ],
    "total_mass_detector": [
        "total_mass", "total_mass_detector", "mtotal", "detector_total_mass",
    ],

    # ratios
    "mass_ratio": [
        "mass_ratio", "q", "massratio",
    ],
    "symmetric_mass_ratio": [
        "symmetric_mass_ratio", "eta",
    ],

    # spins
    "chi_eff": [
        "chi_eff", "chieff", "effective_spin",
    ],
    "chi_p": [
        "chi_p", "chip", "precessing_spin",
    ],
    "a_1": [
        "a_1", "spin_1", "a1", "spin1",
    ],
    "a_2": [
        "a_2", "spin_2", "a2", "spin2",
    ],
    "tilt_1": [
        "tilt_1", "tilt1",
    ],
    "tilt_2": [
        "tilt_2", "tilt2",
    ],

    # distance / redshift
    "luminosity_distance": [
        "luminosity_distance", "distance", "d_l", "luminosity_distance_mpc",
    ],
    "redshift": [
        "redshift", "z",
    ],
}

PRIMARY_VARIABLES = [
    "chirp_mass_source",
    "chirp_mass_detector",
    "total_mass_source",
    "total_mass_detector",
    "mass_1_source",
    "mass_2_source",
    "mass_1_detector",
    "mass_2_detector",
    "mass_ratio",
    "symmetric_mass_ratio",
    "chi_eff",
    "chi_p",
    "luminosity_distance",
    "redshift",
]

SCAN_VARIABLES_POSITIVE = [
    "chirp_mass_source",
    "chirp_mass_detector",
    "total_mass_source",
    "total_mass_detector",
    "mass_1_source",
    "mass_2_source",
    "mass_1_detector",
    "mass_2_detector",
    "mass_ratio",
    "symmetric_mass_ratio",
    "luminosity_distance",
    "redshift",
]

Y_VARIABLES = [
    "chi_eff",
    "chi_p",
    "mass_ratio",
    "symmetric_mass_ratio",
    "log_chirp_mass_source",
    "log_chirp_mass_detector",
    "log_total_mass_source",
    "log_total_mass_detector",
    "log_mass_1_source",
    "log_mass_2_source",
    "log_mass_1_detector",
    "log_mass_2_detector",
]


# =============================================================================
# Utility
# =============================================================================

def safe_name(s: str, max_len: int = 160) -> str:
    out = str(s)
    out = re.sub(r'[<>:"/\\|?*\(\)\[\]\{\}\s]+', "_", out)
    out = re.sub(r"_+", "_", out).strip("._-")
    return (out or "unnamed")[:max_len]


def event_name_from_file(path: str | Path) -> str:
    name = Path(path).name
    m = re.search(r"(GW\d{6}_\d{6})", name)
    if m:
        return m.group(1)
    return Path(path).stem


def md5_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(chunk_size), b""):
            h.update(block)
    return h.hexdigest()


def parse_manifest(path: str | Path) -> pd.DataFrame:
    rows = []
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)

    for line in p.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(maxsplit=1)
        if len(parts) != 2:
            continue
        md5, fname = parts
        rows.append({"md5": md5, "filename": fname})

    return pd.DataFrame(rows)


def verify_manifest(manifest_path: str | Path, base_dir: str | Path = ".") -> pd.DataFrame:
    manifest = parse_manifest(manifest_path)
    base = Path(base_dir)

    rows = []
    for _, r in manifest.iterrows():
        fname = r["filename"]
        expected = r["md5"]
        p = base / fname
        exists = p.exists()
        got = ""
        ok = False
        size_bytes = np.nan
        if exists:
            size_bytes = p.stat().st_size
            got = md5_file(p)
            ok = got.lower() == expected.lower()
        rows.append({
            "filename": fname,
            "expected_md5": expected,
            "exists": bool(exists),
            "actual_md5": got,
            "ok": bool(ok),
            "size_bytes": size_bytes,
        })

    return pd.DataFrame(rows)


def rational_approx(x: float, max_den: int = 64) -> Tuple[str, float, float]:
    if not np.isfinite(x):
        return "", np.nan, np.nan
    f = Fraction(float(x)).limit_denominator(max_den)
    val = f.numerator / f.denominator
    return f"{f.numerator}/{f.denominator}", float(val), float(abs(x - val))


def weighted_lstsq(X: np.ndarray, y: np.ndarray, w: Optional[np.ndarray] = None):
    if w is None:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        chi2 = float(np.sum(resid * resid))
        return beta, chi2, resid

    w = np.asarray(w, dtype=float)
    sw = np.sqrt(np.maximum(w, EPS))
    Xw = X * sw[:, None]
    yw = y * sw
    beta, *_ = np.linalg.lstsq(Xw, yw, rcond=None)
    resid = y - X @ beta
    chi2 = float(np.sum(w * resid * resid))
    return beta, chi2, resid


def make_poly_design(x: np.ndarray, degree: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    z = x - np.nanmean(x)
    scale = np.nanstd(z)
    if np.isfinite(scale) and scale > 0:
        z = z / scale
    return np.column_stack([z ** d for d in range(degree + 1)])


def p_value_ge(real: float, null_vals: np.ndarray) -> float:
    return float((1 + np.sum(null_vals >= real)) / (1 + len(null_vals)))


def tail_count_ge(real: float, null_vals: np.ndarray) -> int:
    return int(np.sum(null_vals >= real))


def robust_width_from_quantiles(samples: np.ndarray) -> float:
    samples = np.asarray(samples, dtype=float)
    samples = samples[np.isfinite(samples)]
    if len(samples) < 5:
        return np.nan
    q16, q84 = np.quantile(samples, [0.16, 0.84])
    return float(0.5 * (q84 - q16))


# =============================================================================
# HDF5 inspection / extraction
# =============================================================================

def hdf5_inventory(path: str | Path, max_items: Optional[int] = None) -> pd.DataFrame:
    rows = []
    path = Path(path)

    with h5py.File(path, "r") as f:
        def visitor(name, obj):
            if max_items is not None and len(rows) >= max_items:
                return
            kind = "group" if isinstance(obj, h5py.Group) else "dataset" if isinstance(obj, h5py.Dataset) else type(obj).__name__
            shape = getattr(obj, "shape", "")
            dtype = str(getattr(obj, "dtype", ""))
            rows.append({
                "file": str(path),
                "name": name,
                "kind": kind,
                "shape": str(shape),
                "dtype": dtype,
            })
        f.visititems(visitor)

    return pd.DataFrame(rows)


def collect_datasets(path: str | Path) -> Dict[str, Dict[str, Any]]:
    """
    Return map:
        full_hdf5_path -> {basename, shape, dtype}
    for numeric datasets.
    """
    out = {}
    path = Path(path)

    with h5py.File(path, "r") as f:
        def visitor(name, obj):
            if not isinstance(obj, h5py.Dataset):
                return

            shape = getattr(obj, "shape", None)
            dtype = getattr(obj, "dtype", None)

            if shape is None or dtype is None:
                return
            if len(shape) == 0:
                return

            # Only numeric arrays.
            try:
                if not np.issubdtype(dtype, np.number):
                    return
            except Exception:
                return

            out[name] = {
                "basename": name.split("/")[-1],
                "shape": tuple(shape),
                "dtype": str(dtype),
            }

        f.visititems(visitor)

    return out


def find_dataset_path(datasets: Dict[str, Dict[str, Any]], aliases: List[str]) -> Optional[str]:
    """
    Heuristic matching:
      1. exact basename match
      2. lowercase exact basename
      3. suffix match
      4. contains match
    """
    aliases_l = [a.lower() for a in aliases]

    # Exact basename.
    for full, info in datasets.items():
        base_l = info["basename"].lower()
        if base_l in aliases_l:
            return full

    # Full path exact/suffix.
    for full, info in datasets.items():
        full_l = full.lower()
        for a in aliases_l:
            if full_l.endswith("/" + a) or full_l == a:
                return full

    # Contains, but avoid choosing covariance or prior if possible.
    preferred = []
    fallback = []
    for full, info in datasets.items():
        full_l = full.lower()
        base_l = info["basename"].lower()
        for a in aliases_l:
            if a in base_l or a in full_l:
                if any(bad in full_l for bad in ["prior", "injection", "cov", "psd"]):
                    fallback.append(full)
                else:
                    preferred.append(full)

    if preferred:
        # Prefer 1D large arrays.
        preferred = sorted(preferred, key=lambda p: (len(datasets[p]["shape"]) != 1, len(p)))
        return preferred[0]
    if fallback:
        return sorted(fallback, key=lambda p: (len(datasets[p]["shape"]) != 1, len(p)))[0]

    return None


def read_numeric_dataset(file_path: str | Path, dataset_path: str, max_samples: Optional[int] = None) -> np.ndarray:
    with h5py.File(file_path, "r") as f:
        arr = np.asarray(f[dataset_path])

    arr = np.asarray(arr)
    arr = np.ravel(arr)
    if max_samples is not None and len(arr) > max_samples:
        # deterministic subsample
        idx = np.linspace(0, len(arr) - 1, max_samples).astype(int)
        arr = arr[idx]
    return arr.astype(float, copy=False)


def derived_from_m1_m2(m1: np.ndarray, m2: np.ndarray) -> Dict[str, np.ndarray]:
    out = {}
    m1 = np.asarray(m1, dtype=float)
    m2 = np.asarray(m2, dtype=float)
    good = np.isfinite(m1) & np.isfinite(m2) & (m1 > 0) & (m2 > 0)
    if not np.any(good):
        return out

    # Ensure m1 >= m2 per sample.
    a = np.maximum(m1, m2)
    b = np.minimum(m1, m2)
    mtot = a + b
    eta = (a * b) / np.maximum(mtot * mtot, EPS)
    q = b / np.maximum(a, EPS)
    mc = (a * b) ** (3.0 / 5.0) / np.maximum(mtot, EPS) ** (1.0 / 5.0)

    out["total_mass_derived"] = mtot
    out["mass_ratio_derived"] = q
    out["symmetric_mass_ratio_derived"] = eta
    out["chirp_mass_derived"] = mc
    return out


def summarize_samples(samples: np.ndarray) -> Dict[str, float]:
    samples = np.asarray(samples, dtype=float)
    samples = samples[np.isfinite(samples)]

    if len(samples) == 0:
        return {
            "median": np.nan, "mean": np.nan, "std": np.nan,
            "q05": np.nan, "q16": np.nan, "q84": np.nan, "q95": np.nan,
            "width_16_84": np.nan, "n_samples": 0,
        }

    q05, q16, q50, q84, q95 = np.quantile(samples, [0.05, 0.16, 0.50, 0.84, 0.95])
    return {
        "median": float(q50),
        "mean": float(np.mean(samples)),
        "std": float(np.std(samples)),
        "q05": float(q05),
        "q16": float(q16),
        "q84": float(q84),
        "q95": float(q95),
        "width_16_84": float(0.5 * (q84 - q16)),
        "n_samples": int(len(samples)),
    }


def extract_event_summary(file_path: str | Path, max_samples: Optional[int] = None) -> Dict[str, Any]:
    """
    Extract medians/widths for known parameters from one event HDF5.
    """
    file_path = Path(file_path)
    event = event_name_from_file(file_path)
    datasets = collect_datasets(file_path)

    row: Dict[str, Any] = {
        "event": event,
        "file": str(file_path),
        "n_numeric_datasets": len(datasets),
    }

    found_arrays: Dict[str, np.ndarray] = {}
    found_paths: Dict[str, str] = {}

    # Direct aliases.
    for canonical, aliases in PARAM_ALIASES.items():
        ds_path = find_dataset_path(datasets, aliases)
        if ds_path is None:
            continue
        try:
            arr = read_numeric_dataset(file_path, ds_path, max_samples=max_samples)
        except Exception as exc:
            row[f"{canonical}_read_error"] = str(exc)
            continue

        arr = arr[np.isfinite(arr)]
        if len(arr) == 0:
            continue

        found_arrays[canonical] = arr
        found_paths[canonical] = ds_path

        s = summarize_samples(arr)
        for key, val in s.items():
            row[f"{canonical}_{key}"] = val
        row[f"{canonical}_path"] = ds_path

    # Derived masses from source if available, else detector.
    for prefix in ["source", "detector"]:
        k1 = f"mass_1_{prefix}"
        k2 = f"mass_2_{prefix}"
        if k1 in found_arrays and k2 in found_arrays:
            derived = derived_from_m1_m2(found_arrays[k1], found_arrays[k2])
            for dname, arr in derived.items():
                canonical = dname.replace("_derived", f"_{prefix}_derived")
                s = summarize_samples(arr)
                for key, val in s.items():
                    row[f"{canonical}_{key}"] = val
            break

    # Fallback canonical variables from derived.
    # Fill missing chirp/total/q/eta from source/detector derived if direct is missing.
    for prefix in ["source", "detector"]:
        for var_base, derived_name in [
            ("chirp_mass", f"chirp_mass_{prefix}_derived"),
            ("total_mass", f"total_mass_{prefix}_derived"),
            ("mass_ratio", f"mass_ratio_{prefix}_derived"),
            ("symmetric_mass_ratio", f"symmetric_mass_ratio_{prefix}_derived"),
        ]:
            target = f"{var_base}_{prefix}" if var_base in ["chirp_mass", "total_mass"] else var_base
            direct_col = f"{target}_median"
            derived_col = f"{derived_name}_median"
            if direct_col not in row and derived_col in row:
                for stat in ["median", "mean", "std", "q05", "q16", "q84", "q95", "width_16_84", "n_samples"]:
                    row[f"{target}_{stat}"] = row.get(f"{derived_name}_{stat}", np.nan)
                row[f"{target}_path"] = f"derived_from_mass_1_{prefix},mass_2_{prefix}"

    return row


# =============================================================================
# Catalog construction
# =============================================================================

def add_log_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in list(df.columns):
        if not col.endswith("_median"):
            continue
        base = col[:-len("_median")]
        vals = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
        if np.nanmax(vals) > 0:
            log_col = f"log_{base}_median"
            df[log_col] = np.where(vals > 0, np.log(vals), np.nan)
    return df


def resolve_variable_column(df: pd.DataFrame, var: str) -> Optional[str]:
    candidates = [
        f"{var}_median",
        var,
    ]

    # if user says log_x, map to log_x_median
    if var.startswith("log_"):
        candidates.insert(0, f"{var}_median")

    for c in candidates:
        if c in df.columns:
            return c

    # handle log_chirp_mass_source -> log_chirp_mass_source_median
    c = f"{var}_median"
    if c in df.columns:
        return c

    return None


def make_scan_dataset(
    event_df: pd.DataFrame,
    x_var: str,
    y_var: str,
    order_by: str = "x",
    min_events: int = 12,
) -> Optional[pd.DataFrame]:
    x_col = resolve_variable_column(event_df, x_var)
    y_col = resolve_variable_column(event_df, y_var)

    if x_col is None or y_col is None:
        return None

    x = pd.to_numeric(event_df[x_col], errors="coerce")
    y = pd.to_numeric(event_df[y_col], errors="coerce")

    # Use y width as optional weight if available.
    y_base = y_col[:-len("_median")] if y_col.endswith("_median") else y_col
    w_col = f"{y_base}_width_16_84"
    if w_col in event_df.columns:
        width = pd.to_numeric(event_df[w_col], errors="coerce").to_numpy(float)
        w = 1.0 / np.maximum(width * width, EPS)
    else:
        w = np.ones(len(event_df), dtype=float)

    out = pd.DataFrame({
        "event": event_df["event"].astype(str),
        "x": x,
        "y": y,
        "w": w,
        "x_var": x_var,
        "y_var": y_var,
        "x_col": x_col,
        "y_col": y_col,
    })

    out = out[np.isfinite(out["x"]) & np.isfinite(out["y"]) & np.isfinite(out["w"]) & (out["w"] > 0)].copy()

    # x must have spread.
    if len(out) < min_events:
        return None
    if out["x"].std() <= 0 or out["y"].std() <= 0:
        return None

    if order_by == "x":
        out = out.sort_values("x").reset_index(drop=True)
    elif order_by == "event":
        out = out.sort_values("event").reset_index(drop=True)
    elif order_by == "y":
        out = out.sort_values("y").reset_index(drop=True)
    else:
        out = out.sort_values("x").reset_index(drop=True)

    return out


# =============================================================================
# Harmonic scan
# =============================================================================

def harmonic_scan(
    x: np.ndarray,
    y: np.ndarray,
    w: Optional[np.ndarray],
    k_grid: np.ndarray,
    degree: int,
) -> Tuple[pd.DataFrame, Dict[str, float], np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if w is not None:
        w = np.asarray(w, dtype=float)

    X0 = make_poly_design(x, degree)
    beta0, chi0, resid0 = weighted_lstsq(X0, y, w=w)

    rows = []
    best = None

    for k in k_grid:
        c = np.cos(k * x)
        s = np.sin(k * x)
        X = np.column_stack([X0, c, s])
        beta, chi, resid = weighted_lstsq(X, y, w=w)
        delta = chi0 - chi
        amp = math.sqrt(beta[-2] ** 2 + beta[-1] ** 2)
        phase = math.atan2(-beta[-1], beta[-2])
        rows.append((k, delta, chi, chi0, amp, phase))

        if best is None or delta > best["delta_chi2"]:
            best = {
                "k_best": float(k),
                "delta_chi2": float(delta),
                "chi2_baseline": float(chi0),
                "chi2_harmonic": float(chi),
                "amplitude": float(amp),
                "phase": float(phase),
                "baseline_degree": int(degree),
                "n_points": int(len(x)),
            }

    scan_df = pd.DataFrame(rows, columns=[
        "k", "delta_chi2", "chi2_harmonic", "chi2_baseline", "amplitude", "phase"
    ])

    return scan_df, best, resid0, X0 @ beta0


def null_scan(
    x: np.ndarray,
    y: np.ndarray,
    w: Optional[np.ndarray],
    k_grid: np.ndarray,
    degree: int,
    real_delta: float,
    null_n: int,
    seed: int,
) -> Tuple[np.ndarray, float, int]:
    rng = np.random.default_rng(seed)

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    w = None if w is None else np.asarray(w, dtype=float)

    X0 = make_poly_design(x, degree)
    beta0, chi0, resid0 = weighted_lstsq(X0, y, w=w)
    baseline = X0 @ beta0

    null_max = np.empty(null_n, dtype=float)

    for i in range(null_n):
        y_null = baseline + resid0[rng.permutation(len(resid0))]
        scan_df, best, _, _ = harmonic_scan(x, y_null, w, k_grid, degree)
        null_max[i] = best["delta_chi2"]

        if (i + 1) % 500 == 0 or i + 1 == null_n:
            print(f"[null] {i+1}/{null_n}")

    p = p_value_ge(real_delta, null_max)
    count = tail_count_ge(real_delta, null_max)
    return null_max, p, count


# =============================================================================
# Branch / ratio matching
# =============================================================================

def delta_domain(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.nanmax(x) - np.nanmin(x))


def n_from_k(k: float, delta: float) -> float:
    return float(k * delta / (2.0 * math.pi))


def k_from_n(n: float, delta: float) -> float:
    return float(2.0 * math.pi * n / delta)


def branch_grid_matches(k_best: float, delta: float, max_n: int = 40) -> pd.DataFrame:
    n_obs = n_from_k(k_best, delta)
    rows = []

    for n in range(1, max_n + 1):
        k_n = k_from_n(n, delta)
        rows.append({
            "n": n,
            "k_n": k_n,
            "k_best": k_best,
            "n_obs": n_obs,
            "k_error": k_best - k_n,
            "abs_k_error": abs(k_best - k_n),
            "n_error": n_obs - n,
            "abs_n_error": abs(n_obs - n),
        })

    return pd.DataFrame(rows).sort_values("abs_k_error").reset_index(drop=True)


def ratio_scan_results(scan_df: pd.DataFrame, delta: float, top_peaks: int = 12) -> pd.DataFrame:
    """
    Find prominent local maxima in scan curve and compare peak ratios to branch ratios.
    """
    df = scan_df.copy()
    if df.empty:
        return pd.DataFrame()

    y = df["delta_chi2"].to_numpy(float)
    if find_peaks is not None and len(y) >= 5:
        peaks, props = find_peaks(y, distance=max(1, len(y) // 80))
        if len(peaks) == 0:
            peaks = np.array([int(np.argmax(y))])
    else:
        peaks = np.array([int(np.argmax(y))])

    peak_df = df.iloc[peaks].copy()
    peak_df = peak_df.sort_values("delta_chi2", ascending=False).head(top_peaks).reset_index(drop=True)

    rows = []
    for i in range(len(peak_df)):
        for j in range(i + 1, len(peak_df)):
            k1 = float(peak_df.loc[i, "k"])
            k2 = float(peak_df.loc[j, "k"])
            lo, hi = min(k1, k2), max(k1, k2)
            ratio = hi / lo
            rat, rat_val, rat_err = rational_approx(ratio, max_den=64)

            # Compare to named branch ratios.
            branch_best = None
            for bname, ns in BRANCHES.items():
                ns_sorted = sorted(set(float(n) for n in ns if n > 0))
                for a in range(len(ns_sorted)):
                    for b in range(a + 1, len(ns_sorted)):
                        br = ns_sorted[b] / ns_sorted[a]
                        err = abs(ratio - br)
                        if branch_best is None or err < branch_best["branch_ratio_error"]:
                            branch_best = {
                                "nearest_branch_name": bname,
                                "nearest_branch_low": ns_sorted[a],
                                "nearest_branch_high": ns_sorted[b],
                                "nearest_branch_ratio": br,
                                "branch_ratio_error": err,
                            }

            row = {
                "k_low": lo,
                "k_high": hi,
                "ratio_high_low": ratio,
                "ratio_rational": rat,
                "ratio_rational_error": rat_err,
                "delta_low": float(peak_df.loc[j if k2 == lo else i, "delta_chi2"]),
                "delta_high": float(peak_df.loc[i if k1 == hi else j, "delta_chi2"]),
            }
            if branch_best:
                row.update(branch_best)
            rows.append(row)

    if not rows:
        return pd.DataFrame()

    return pd.DataFrame(rows).sort_values(
        ["branch_ratio_error", "ratio_rational_error"],
        ascending=[True, True],
    ).reset_index(drop=True)


# =============================================================================
# Full pipeline
# =============================================================================

def discover_files(input_patterns: List[str]) -> List[str]:
    files = []
    for pat in input_patterns:
        hits = glob.glob(pat)
        if hits:
            files.extend(hits)
        elif Path(pat).exists():
            files.append(pat)
    return sorted(set(files))


def build_event_table(files: List[str], max_samples: Optional[int], out_dir: Path) -> pd.DataFrame:
    rows = []
    errors = []

    for i, f in enumerate(files, 1):
        print(f"[extract] {i}/{len(files)} {f}")
        try:
            row = extract_event_summary(f, max_samples=max_samples)
            rows.append(row)
        except Exception as exc:
            print(f"[warn] failed {f}: {exc}")
            errors.append({"file": f, "error": str(exc)})

    df = pd.DataFrame(rows)
    if not df.empty:
        df = add_log_columns(df)
        df.to_csv(out_dir / "gwtc4_event_summary.csv", index=False)

    if errors:
        pd.DataFrame(errors).to_csv(out_dir / "gwtc4_extract_errors.csv", index=False)

    return df


def choose_variable_pairs(event_df: pd.DataFrame, requested_pairs: Optional[List[str]] = None) -> List[Tuple[str, str]]:
    if requested_pairs:
        pairs = []
        for item in requested_pairs:
            if ":" not in item:
                raise ValueError(f"Pair must be X:Y, got {item}")
            a, b = item.split(":", 1)
            pairs.append((a.strip(), b.strip()))
        return pairs

    # Conservative default pairs:
    # x is a positive log mass/distance/ration coordinate;
    # y is a shape/spin/other coordinate.
    pairs = []
    x_candidates = []
    for v in SCAN_VARIABLES_POSITIVE:
        logv = f"log_{v}"
        if resolve_variable_column(event_df, logv) is not None:
            x_candidates.append(logv)
        elif resolve_variable_column(event_df, v) is not None:
            x_candidates.append(v)

    y_candidates = []
    for v in Y_VARIABLES:
        if resolve_variable_column(event_df, v) is not None:
            y_candidates.append(v)

    # Keep it manageable and meaningful.
    for x in x_candidates:
        for y in y_candidates:
            if x == y:
                continue
            pairs.append((x, y))

    # Also scan index/order against log variables as y; use event order x_index.
    # This is added later separately if needed.
    return pairs


def run_catalog_scans(
    event_df: pd.DataFrame,
    out_dir: Path,
    variable_pairs: List[Tuple[str, str]],
    k_grid: np.ndarray,
    degree: int,
    null_n: int,
    min_events: int,
    seed: int,
    skip_null: bool = False,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    scans_dir = out_dir / "scans"
    nulls_dir = out_dir / "nulls"
    scans_dir.mkdir(exist_ok=True, parents=True)
    nulls_dir.mkdir(exist_ok=True, parents=True)

    scan_rows = []
    branch_rows = []
    ratio_rows = []

    rng = np.random.default_rng(seed)

    for idx, (x_var, y_var) in enumerate(variable_pairs):
        ds = make_scan_dataset(event_df, x_var, y_var, min_events=min_events)
        if ds is None:
            continue

        x = ds["x"].to_numpy(float)
        y = ds["y"].to_numpy(float)
        w = ds["w"].to_numpy(float)

        if len(x) < min_events:
            continue

        label = safe_name(f"{x_var}__{y_var}")
        print("\n" + "=" * 80)
        print(f"[scan] {x_var} -> {y_var}  N={len(x)}")
        print("=" * 80)

        scan_df, best, resid0, baseline = harmonic_scan(x, y, w, k_grid, degree)
        scan_path = scans_dir / f"scan_{label}.csv"
        scan_df.to_csv(scan_path, index=False)

        ddom = delta_domain(x)
        n_obs = n_from_k(best["k_best"], ddom)

        best.update({
            "x_var": x_var,
            "y_var": y_var,
            "n_events": int(len(x)),
            "x_min": float(np.min(x)),
            "x_max": float(np.max(x)),
            "delta_x": float(ddom),
            "n_obs": float(n_obs),
            "scan_csv": str(scan_path),
        })

        if not skip_null and null_n > 0:
            null_seed = int(rng.integers(0, 2**31 - 1))
            null_max, p, count = null_scan(
                x=x,
                y=y,
                w=w,
                k_grid=k_grid,
                degree=degree,
                real_delta=best["delta_chi2"],
                null_n=null_n,
                seed=null_seed,
            )
            null_path = nulls_dir / f"null_{label}.csv"
            pd.DataFrame({"null_max_delta_chi2": null_max}).to_csv(null_path, index=False)

            best.update({
                "scan_null_p": float(p),
                "tail_count_ge": int(count),
                "null_n": int(null_n),
                "null_max_median": float(np.median(null_max)),
                "null_max_p95": float(np.quantile(null_max, 0.95)),
                "null_max_p99": float(np.quantile(null_max, 0.99)),
                "null_csv": str(null_path),
            })
        else:
            best.update({
                "scan_null_p": np.nan,
                "tail_count_ge": np.nan,
                "null_n": 0,
                "null_max_median": np.nan,
                "null_max_p95": np.nan,
                "null_max_p99": np.nan,
                "null_csv": "",
            })

        scan_rows.append(best)

        # Branch nearest integer grid.
        bm = branch_grid_matches(best["k_best"], ddom, max_n=40).head(10)
        bm["x_var"] = x_var
        bm["y_var"] = y_var
        bm["scan_k_best"] = best["k_best"]
        bm["scan_delta_chi2"] = best["delta_chi2"]
        bm["scan_null_p"] = best["scan_null_p"]
        branch_rows.append(bm)

        # Peak ratio candidates.
        rr = ratio_scan_results(scan_df, ddom, top_peaks=12)
        if not rr.empty:
            rr["x_var"] = x_var
            rr["y_var"] = y_var
            rr["scan_k_best"] = best["k_best"]
            rr["scan_delta_chi2"] = best["delta_chi2"]
            rr["scan_null_p"] = best["scan_null_p"]
            ratio_rows.append(rr)

        print(json.dumps({k: best[k] for k in [
            "x_var", "y_var", "n_events", "k_best", "delta_chi2",
            "n_obs", "scan_null_p", "tail_count_ge", "null_n"
        ]}, indent=2))

    scans = pd.DataFrame(scan_rows)
    branches = pd.concat(branch_rows, ignore_index=True) if branch_rows else pd.DataFrame()
    ratios = pd.concat(ratio_rows, ignore_index=True) if ratio_rows else pd.DataFrame()

    if not scans.empty:
        scans = scans.sort_values(
            ["scan_null_p", "delta_chi2"],
            ascending=[True, False],
            na_position="last",
        ).reset_index(drop=True)
        scans.to_csv(out_dir / "gwtc4_scan_results.csv", index=False)

    if not branches.empty:
        branches = branches.sort_values(
            ["scan_null_p", "abs_k_error"],
            ascending=[True, True],
            na_position="last",
        ).reset_index(drop=True)
        branches.to_csv(out_dir / "gwtc4_branch_match_results.csv", index=False)

    if not ratios.empty:
        ratios = ratios.sort_values(
            ["scan_null_p", "branch_ratio_error", "ratio_rational_error"],
            ascending=[True, True, True],
            na_position="last",
        ).reset_index(drop=True)
        ratios.to_csv(out_dir / "gwtc4_ratio_results.csv", index=False)

    return scans, branches, ratios


def inspect_files(files: List[str], out_dir: Path, limit: int = 3, max_items: int = 300):
    rows = []
    for f in files[:limit]:
        print(f"[inspect] {f}")
        try:
            inv = hdf5_inventory(f, max_items=max_items)
            rows.append(inv)
            print(inv.head(80).to_string(index=False))
        except Exception as exc:
            print(f"[warn] inspect failed {f}: {exc}")

    if rows:
        out = pd.concat(rows, ignore_index=True)
        out.to_csv(out_dir / "hdf5_inventory.csv", index=False)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", nargs="*", default=["*combined_PEDataRelease.hdf5"],
                        help="Input HDF5 glob(s). Default: *combined_PEDataRelease.hdf5")
    parser.add_argument("--manifest", default=None, help="Optional md5sums.txt")
    parser.add_argument("--base-dir", default=".", help="Base dir for manifest verification")
    parser.add_argument("--out-dir", default=OUTDIR_DEFAULT)

    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--inspect-only", action="store_true")
    parser.add_argument("--inspect-limit", type=int, default=3)
    parser.add_argument("--inspect-max-items", type=int, default=300)

    parser.add_argument("--max-samples", type=int, default=None,
                        help="Optional deterministic subsample per posterior parameter.")
    parser.add_argument("--min-events", type=int, default=MIN_EVENTS_DEFAULT)

    parser.add_argument("--k-min", type=float, default=K_MIN_DEFAULT)
    parser.add_argument("--k-max", type=float, default=K_MAX_DEFAULT)
    parser.add_argument("--n-k", type=int, default=N_K_DEFAULT)
    parser.add_argument("--degree", type=int, default=BASELINE_DEGREE_DEFAULT)
    parser.add_argument("--null-n", type=int, default=NULL_N_DEFAULT)
    parser.add_argument("--skip-null", action="store_true")
    parser.add_argument("--seed", type=int, default=RNG_SEED_DEFAULT)

    parser.add_argument("--pair", action="append", default=None,
                        help="Specific scan pair X:Y. Can repeat. Example: --pair log_chirp_mass_source:chi_eff")

    parser.add_argument("--top-n", type=int, default=30)

    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    run_summary: Dict[str, Any] = {
        "args": vars(args),
        "delta_ell_lhc": DELTA_ELL_LHC,
        "branches": BRANCHES,
    }

    # Verify manifest.
    if args.manifest:
        print(f"[manifest] verifying {args.manifest}")
        vf = verify_manifest(args.manifest, base_dir=args.base_dir)
        vf.to_csv(out_dir / "manifest_verification.csv", index=False)
        run_summary["manifest_total"] = int(len(vf))
        run_summary["manifest_present"] = int(vf["exists"].sum())
        run_summary["manifest_ok"] = int(vf["ok"].sum())
        print(vf["ok"].value_counts(dropna=False).to_string())
        print(f"[save] {out_dir / 'manifest_verification.csv'}")

        if args.verify_only:
            with open(out_dir / "gwtc4_run_summary.json", "w", encoding="utf-8") as f:
                json.dump(run_summary, f, indent=2)
            return

    files = discover_files(args.input)
    files = [f for f in files if f.endswith(".hdf5") or f.endswith(".h5")]
    files = sorted(files)

    if not files:
        raise FileNotFoundError(f"No HDF5 files found for input: {args.input}")

    # Exclude PESummaryTable from per-event extraction unless user explicitly only passed it.
    event_files = [f for f in files if "PESummaryTable" not in Path(f).name]
    if not event_files:
        event_files = files

    print(f"[files] {len(event_files)} event HDF5 files")
    for f in event_files[:20]:
        print("  ", f)
    if len(event_files) > 20:
        print(f"  ... {len(event_files)-20} more")

    run_summary["n_files"] = len(event_files)
    run_summary["files"] = event_files

    if args.inspect_only:
        inspect_files(event_files, out_dir, limit=args.inspect_limit, max_items=args.inspect_max_items)
        with open(out_dir / "gwtc4_run_summary.json", "w", encoding="utf-8") as f:
            json.dump(run_summary, f, indent=2)
        return

    # Extract event table.
    event_df = build_event_table(event_files, max_samples=args.max_samples, out_dir=out_dir)
    if event_df.empty:
        raise RuntimeError("No event summaries could be extracted.")

    print(f"[events] extracted {len(event_df)} event rows")
    run_summary["n_events_extracted"] = int(len(event_df))
    run_summary["event_columns"] = list(event_df.columns)

    # Build variable pairs.
    variable_pairs = choose_variable_pairs(event_df, requested_pairs=args.pair)
    if not variable_pairs:
        raise RuntimeError("No valid variable pairs found. Try --inspect-only to see HDF5 paths.")

    print(f"[pairs] {len(variable_pairs)} candidate scan pairs")

    k_grid = np.linspace(args.k_min, args.k_max, args.n_k)

    scans, branches, ratios = run_catalog_scans(
        event_df=event_df,
        out_dir=out_dir,
        variable_pairs=variable_pairs,
        k_grid=k_grid,
        degree=args.degree,
        null_n=args.null_n,
        min_events=args.min_events,
        seed=args.seed,
        skip_null=args.skip_null,
    )

    run_summary["n_scans"] = int(len(scans))
    run_summary["n_branch_rows"] = int(len(branches))
    run_summary["n_ratio_rows"] = int(len(ratios))

    if not scans.empty:
        run_summary["top_scans"] = scans.head(args.top_n).to_dict(orient="records")
    if not branches.empty:
        run_summary["top_branch_matches"] = branches.head(args.top_n).to_dict(orient="records")
    if not ratios.empty:
        run_summary["top_ratio_matches"] = ratios.head(args.top_n).to_dict(orient="records")

    with open(out_dir / "gwtc4_run_summary.json", "w", encoding="utf-8") as f:
        json.dump(run_summary, f, indent=2)

    print("\n" + "#" * 100)
    print("[done]")
    print("#" * 100)
    print(f"[save] {out_dir / 'gwtc4_event_summary.csv'}")
    print(f"[save] {out_dir / 'gwtc4_scan_results.csv'}")
    print(f"[save] {out_dir / 'gwtc4_branch_match_results.csv'}")
    print(f"[save] {out_dir / 'gwtc4_ratio_results.csv'}")
    print(f"[save] {out_dir / 'gwtc4_run_summary.json'}")

    if not scans.empty:
        print("\n[top scans]")
        cols = ["x_var", "y_var", "n_events", "k_best", "delta_chi2", "n_obs", "scan_null_p", "tail_count_ge", "null_n"]
        print(scans.head(args.top_n)[cols].to_string(index=False))

    if not branches.empty:
        print("\n[top branch nearest-integer matches]")
        cols = ["x_var", "y_var", "n", "k_n", "k_best", "abs_k_error", "n_obs", "abs_n_error", "scan_null_p"]
        print(branches.head(args.top_n)[cols].to_string(index=False))

    if not ratios.empty:
        print("\n[top peak-ratio matches]")
        cols = ["x_var", "y_var", "k_low", "k_high", "ratio_high_low", "ratio_rational",
                "nearest_branch_name", "nearest_branch_low", "nearest_branch_high",
                "nearest_branch_ratio", "branch_ratio_error", "scan_null_p"]
        print(ratios.head(args.top_n)[cols].to_string(index=False))


if __name__ == "__main__":
    main()
