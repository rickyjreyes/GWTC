#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
GWTC-4 WCT Posterior Scanner — Compound posterior_samples version
-----------------------------------------------------------------

This version fixes the issue found by the pathfinder output:

    GWTC-4 posterior samples are stored as COMPOUND HDF5 datasets like
        C00:Mixed+XO4a/posterior_samples
        C00:Mixed/posterior_samples
        C00:SEOBNRv5PHM/posterior_samples

Each dataset is a table-like array whose dtype fields contain:
    chirp_mass_source, mass_1_source, mass_2_source, mass_ratio,
    symmetric_mass_ratio, chi_eff, chi_p, luminosity_distance, redshift, ...

The earlier posterior-only scanner looked for one parameter per HDF5 dataset,
so it missed these compound posterior tables. This script reads fields from
the compound posterior_samples table directly.

Install:
    pip install numpy pandas scipy h5py

Fast test:
    python gwtc4_wct_catalog_scan_compound_posterior.py --input "*combined_PEDataRelease.hdf5" --null-n 500 --pair log_chirp_mass_source:symmetric_mass_ratio

Full focused test:
    python gwtc4_wct_catalog_scan_compound_posterior.py --input "*combined_PEDataRelease.hdf5" --null-n 5000 --k-max 80 --pair log_chirp_mass_source:symmetric_mass_ratio

Outputs:
    outputs_gwtc4_wct_compound/
        gwtc4_event_summary.csv
        gwtc4_posterior_table_audit.csv
        gwtc4_scan_results.csv
        gwtc4_branch_match_results.csv
        gwtc4_ratio_results.csv
        gwtc4_run_summary.json
        scans/
        nulls/
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import re
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


EPS = 1e-12
OUTDIR_DEFAULT = "outputs_gwtc4_wct_compound"

K_MIN_DEFAULT = 0.5
K_MAX_DEFAULT = 80.0
N_K_DEFAULT = 4000
NULL_N_DEFAULT = 5000
BASELINE_DEGREE_DEFAULT = 2
MIN_EVENTS_DEFAULT = 12
SEED_DEFAULT = 12345

DELTA_ELL_LHC = 4.780150335923678

POSTERIOR_DATASET_NAME = "posterior_samples"

PREFERRED_POSTERIOR_PATH_SUBSTRINGS = [
    "Mixed+XO4a/posterior_samples",
    "Mixed/posterior_samples",
    "Mixed:HighSpin/posterior_samples",
    "Mixed:LowSpinSecondary/posterior_samples",
    "SEOBNRv5PHM/posterior_samples",
    "IMRPhenomXPHM/posterior_samples",
    "IMRPhenomXO4a/posterior_samples",
    "NRSur7dq4/posterior_samples",
]

PARAM_ALIASES = {
    "mass_1_source": ["mass_1_source", "m1_source"],
    "mass_2_source": ["mass_2_source", "m2_source"],
    "mass_1": ["mass_1", "m1"],
    "mass_2": ["mass_2", "m2"],
    "chirp_mass_source": ["chirp_mass_source", "mc_source"],
    "chirp_mass": ["chirp_mass", "mc"],
    "total_mass_source": ["total_mass_source"],
    "total_mass": ["total_mass"],
    "mass_ratio": ["mass_ratio", "q"],
    "inverted_mass_ratio": ["inverted_mass_ratio"],
    "symmetric_mass_ratio": ["symmetric_mass_ratio", "eta"],
    "chi_eff": ["chi_eff", "chieff"],
    "chi_p": ["chi_p", "chip"],
    "a_1": ["a_1"],
    "a_2": ["a_2"],
    "tilt_1": ["tilt_1"],
    "tilt_2": ["tilt_2"],
    "luminosity_distance": ["luminosity_distance"],
    "redshift": ["redshift"],
    "final_mass_source": ["final_mass_source"],
    "final_spin": ["final_spin"],
    "network_matched_filter_snr": ["network_matched_filter_snr"],
    "network_optimal_snr": ["network_optimal_snr"],
    "log_likelihood": ["log_likelihood"],
    "weights": ["weights"],
}

SCAN_VARIABLES_POSITIVE = [
    "chirp_mass_source",
    "chirp_mass",
    "total_mass_source",
    "total_mass",
    "mass_1_source",
    "mass_2_source",
    "mass_1",
    "mass_2",
    "mass_ratio",
    "inverted_mass_ratio",
    "symmetric_mass_ratio",
    "luminosity_distance",
    "redshift",
    "final_mass_source",
]

Y_VARIABLES = [
    "symmetric_mass_ratio",
    "mass_ratio",
    "inverted_mass_ratio",
    "chi_eff",
    "chi_p",
    "a_1",
    "a_2",
    "tilt_1",
    "tilt_2",
    "final_spin",
    "network_matched_filter_snr",
    "network_optimal_snr",
    "log_chirp_mass_source",
    "log_chirp_mass",
    "log_total_mass_source",
    "log_total_mass",
    "log_mass_1_source",
    "log_mass_2_source",
    "log_mass_1",
    "log_mass_2",
    "log_luminosity_distance",
    "log_redshift",
    "log_final_mass_source",
]

BRANCHES = {
    "star_transposed_6_17": [6.0, 17.0],
    "koide_10_15_20": [10.0, 15.0, 20.0],
    "folded_4over9": [20.0 / 3.0, 15.0, 40.0 / 3.0],
    "integer_1_to_40": list(map(float, range(1, 41))),
}


def safe_name(s: str, max_len: int = 180) -> str:
    out = re.sub(r'[<>:"/\\|?*\(\)\[\]\{\}\s]+', "_", str(s))
    out = re.sub(r"_+", "_", out).strip("._-")
    return (out or "unnamed")[:max_len]


def event_name_from_file(path: str | Path) -> str:
    name = Path(path).name
    m = re.search(r"(GW\d{6}_\d{6})", name)
    return m.group(1) if m else Path(path).stem


def discover_files(patterns: List[str]) -> List[str]:
    files = []
    for pat in patterns:
        hits = glob.glob(pat)
        if hits:
            files.extend(hits)
        elif Path(pat).exists():
            files.append(pat)
    return sorted(set(files))


def is_compound_dataset(obj: Any) -> bool:
    return isinstance(obj, h5py.Dataset) and getattr(obj, "dtype", None) is not None and obj.dtype.fields is not None


def posterior_tables(path: str | Path) -> List[Dict[str, Any]]:
    rows = []
    with h5py.File(path, "r") as f:
        def visitor(name: str, obj: Any):
            if not is_compound_dataset(obj):
                return
            if not name.endswith(POSTERIOR_DATASET_NAME):
                return
            low = name.lower()
            if "prior" in low or "priors" in low:
                return
            fields = list(obj.dtype.fields.keys())
            rows.append({
                "file": str(path),
                "event": event_name_from_file(path),
                "dataset": name,
                "n_samples": int(obj.shape[0]) if len(obj.shape) else 0,
                "n_fields": len(fields),
                "fields": ",".join(fields),
                "has_weights": "weights" in fields,
                "score": posterior_table_score(name, fields, int(obj.shape[0]) if len(obj.shape) else 0),
            })
        f.visititems(visitor)
    return sorted(rows, key=lambda r: (-r["score"], -r["n_samples"], r["dataset"]))


def posterior_table_score(name: str, fields: List[str], n_samples: int) -> float:
    score = 0.0
    for i, pref in enumerate(PREFERRED_POSTERIOR_PATH_SUBSTRINGS):
        if pref.lower() in name.lower():
            score += 1000.0 - 10.0 * i
            break
    useful = [
        "chirp_mass_source", "mass_1_source", "mass_2_source", "mass_ratio",
        "symmetric_mass_ratio", "chi_eff", "chi_p", "luminosity_distance", "redshift"
    ]
    score += 5.0 * sum(1 for u in useful if u in fields)
    score += min(float(n_samples) / 1000.0, 100.0)
    return score


def choose_posterior_table(tables: List[Dict[str, Any]], strategy: str = "preferred") -> Optional[Dict[str, Any]]:
    if not tables:
        return None
    if strategy == "largest":
        return sorted(tables, key=lambda r: (-r["n_samples"], r["dataset"]))[0]
    if strategy == "first":
        return sorted(tables, key=lambda r: r["dataset"])[0]
    if strategy == "preferred":
        return sorted(tables, key=lambda r: (-r["score"], -r["n_samples"], r["dataset"]))[0]
    raise ValueError(f"Unknown posterior choice strategy: {strategy}")


def resolve_field(fields: Iterable[str], aliases: List[str]) -> Optional[str]:
    fields = list(fields)
    lower = {f.lower(): f for f in fields}
    for a in aliases:
        if a.lower() in lower:
            return lower[a.lower()]
    for a in aliases:
        al = a.lower()
        for f in fields:
            if al == f.lower() or al in f.lower():
                return f
    return None


def weighted_quantile(values: np.ndarray, quantiles: List[float], weights: Optional[np.ndarray] = None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values)
    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        good &= np.isfinite(weights) & (weights > 0)
    values = values[good]
    if len(values) == 0:
        return np.full(len(quantiles), np.nan)
    if weights is None:
        return np.quantile(values, quantiles)
    weights = weights[good]
    order = np.argsort(values)
    values = values[order]
    weights = weights[order]
    cdf = np.cumsum(weights)
    cdf = cdf / cdf[-1]
    return np.interp(quantiles, cdf, values)


def summarize_samples(values: np.ndarray, weights: Optional[np.ndarray] = None) -> Dict[str, float]:
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values)
    if weights is not None:
        weights = np.asarray(weights, dtype=float)
        good &= np.isfinite(weights) & (weights > 0)
    v = values[good]
    w = weights[good] if weights is not None else None
    if len(v) == 0:
        return {k: np.nan for k in ["median", "mean", "std", "q05", "q16", "q84", "q95", "width_16_84"]} | {"n_samples": 0}
    q05, q16, q50, q84, q95 = weighted_quantile(v, [0.05, 0.16, 0.5, 0.84, 0.95], w)
    if w is None:
        mean = float(np.mean(v))
        std = float(np.std(v))
    else:
        mean = float(np.average(v, weights=w))
        std = float(np.sqrt(np.average((v - mean) ** 2, weights=w)))
    return {
        "median": float(q50),
        "mean": mean,
        "std": std,
        "q05": float(q05),
        "q16": float(q16),
        "q84": float(q84),
        "q95": float(q95),
        "width_16_84": float(0.5 * (q84 - q16)),
        "n_samples": int(len(v)),
    }


def extract_event_summary(file_path: str | Path, posterior_choice: str, max_samples: Optional[int] = None, use_weights: bool = False) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    file_path = Path(file_path)
    tables = posterior_tables(file_path)
    chosen = choose_posterior_table(tables, strategy=posterior_choice)

    row: Dict[str, Any] = {
        "event": event_name_from_file(file_path),
        "file": str(file_path),
        "posterior_dataset": "",
        "posterior_n_samples": 0,
        "posterior_n_fields": 0,
    }

    if chosen is None:
        row["extract_error"] = "no posterior_samples compound dataset found"
        return row, tables

    row["posterior_dataset"] = chosen["dataset"]
    row["posterior_n_samples"] = int(chosen["n_samples"])
    row["posterior_n_fields"] = int(chosen["n_fields"])

    with h5py.File(file_path, "r") as f:
        d = f[chosen["dataset"]]
        fields = list(d.dtype.fields.keys())

        n = d.shape[0]
        if max_samples is not None and n > max_samples:
            idx = np.linspace(0, n - 1, max_samples).astype(int)
            arr = d[idx]
        else:
            arr = d[:]

    weights = None
    if use_weights and "weights" in fields:
        weights = np.asarray(arr["weights"], dtype=float)

    for canonical, aliases in PARAM_ALIASES.items():
        field = resolve_field(fields, aliases)
        if field is None:
            continue
        vals = np.asarray(arr[field], dtype=float)
        s = summarize_samples(vals, weights=weights)
        for k, v in s.items():
            row[f"{canonical}_{k}"] = v
        row[f"{canonical}_path"] = f"{chosen['dataset']}[{field}]"

    return row, tables


def add_log_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in list(df.columns):
        if not col.endswith("_median"):
            continue
        base = col[:-len("_median")]
        vals = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
        if np.any(np.isfinite(vals) & (vals > 0)):
            df[f"log_{base}_median"] = np.where(vals > 0, np.log(vals), np.nan)
    return df


def resolve_variable_column(df: pd.DataFrame, var: str) -> Optional[str]:
    for c in [f"{var}_median", var]:
        if c in df.columns:
            return c
    return None


def make_poly_design(x: np.ndarray, degree: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    z = x - np.nanmean(x)
    s = np.nanstd(z)
    if np.isfinite(s) and s > 0:
        z = z / s
    return np.column_stack([z ** d for d in range(degree + 1)])


def weighted_lstsq(X: np.ndarray, y: np.ndarray, w: Optional[np.ndarray] = None):
    if w is None:
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        return beta, float(np.sum(resid * resid)), resid
    sw = np.sqrt(np.maximum(w, EPS))
    beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    resid = y - X @ beta
    return beta, float(np.sum(w * resid * resid)), resid


def harmonic_scan(x: np.ndarray, y: np.ndarray, w: Optional[np.ndarray], k_grid: np.ndarray, degree: int):
    X0 = make_poly_design(x, degree)
    beta0, chi0, resid0 = weighted_lstsq(X0, y, w)
    rows = []
    best = None
    for k in k_grid:
        X = np.column_stack([X0, np.cos(k * x), np.sin(k * x)])
        beta, chi, resid = weighted_lstsq(X, y, w)
        delta = chi0 - chi
        amp = float(math.sqrt(beta[-2] ** 2 + beta[-1] ** 2))
        phase = float(math.atan2(-beta[-1], beta[-2]))
        rows.append((float(k), float(delta), float(chi), float(chi0), amp, phase))
        if best is None or delta > best["delta_chi2"]:
            best = {
                "k_best": float(k),
                "delta_chi2": float(delta),
                "chi2_baseline": float(chi0),
                "chi2_harmonic": float(chi),
                "amplitude": amp,
                "phase": phase,
                "baseline_degree": int(degree),
                "n_points": int(len(x)),
            }
    scan_df = pd.DataFrame(rows, columns=["k", "delta_chi2", "chi2_harmonic", "chi2_baseline", "amplitude", "phase"])
    return scan_df, best, resid0, X0 @ beta0


def null_scan(x: np.ndarray, y: np.ndarray, w: Optional[np.ndarray], k_grid: np.ndarray, degree: int, real_delta: float, null_n: int, seed: int):
    rng = np.random.default_rng(seed)
    X0 = make_poly_design(x, degree)
    beta0, chi0, resid0 = weighted_lstsq(X0, y, w)
    baseline = X0 @ beta0
    null_max = np.empty(null_n, dtype=float)
    for i in range(null_n):
        y_null = baseline + resid0[rng.permutation(len(resid0))]
        _, best, _, _ = harmonic_scan(x, y_null, w, k_grid, degree)
        null_max[i] = best["delta_chi2"]
        if (i + 1) % 500 == 0 or i + 1 == null_n:
            print(f"[null] {i+1}/{null_n}")
    p = float((1 + np.sum(null_max >= real_delta)) / (1 + len(null_max)))
    count = int(np.sum(null_max >= real_delta))
    return null_max, p, count


def delta_domain(x: np.ndarray) -> float:
    return float(np.nanmax(x) - np.nanmin(x))


def n_from_k(k: float, delta: float) -> float:
    return float(k * delta / (2.0 * math.pi))


def k_from_n(n: float, delta: float) -> float:
    return float(2.0 * math.pi * n / delta)


def branch_grid_matches(k_best: float, delta: float, max_n: int = 60) -> pd.DataFrame:
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


def rational_approx(x: float, max_den: int = 64) -> Tuple[str, float, float]:
    if not np.isfinite(x):
        return "", np.nan, np.nan
    f = Fraction(float(x)).limit_denominator(max_den)
    val = f.numerator / f.denominator
    return f"{f.numerator}/{f.denominator}", float(val), float(abs(x - val))


def ratio_scan_results(scan_df: pd.DataFrame, top_peaks: int = 12) -> pd.DataFrame:
    if scan_df.empty:
        return pd.DataFrame()
    y = scan_df["delta_chi2"].to_numpy(float)
    if find_peaks is not None and len(y) >= 5:
        peaks, _ = find_peaks(y, distance=max(1, len(y) // 80))
        if len(peaks) == 0:
            peaks = np.array([int(np.argmax(y))])
    else:
        peaks = np.array([int(np.argmax(y))])
    peak_df = scan_df.iloc[peaks].sort_values("delta_chi2", ascending=False).head(top_peaks).reset_index(drop=True)

    rows = []
    for i in range(len(peak_df)):
        for j in range(i + 1, len(peak_df)):
            k1 = float(peak_df.loc[i, "k"])
            k2 = float(peak_df.loc[j, "k"])
            lo, hi = min(k1, k2), max(k1, k2)
            ratio = hi / lo
            rat, rat_val, rat_err = rational_approx(ratio)
            branch_best = None
            for bname, ns in BRANCHES.items():
                ns = sorted(set(float(n) for n in ns if n > 0))
                for a in range(len(ns)):
                    for b in range(a + 1, len(ns)):
                        br = ns[b] / ns[a]
                        err = abs(ratio - br)
                        if branch_best is None or err < branch_best["branch_ratio_error"]:
                            branch_best = {
                                "nearest_branch_name": bname,
                                "nearest_branch_low": ns[a],
                                "nearest_branch_high": ns[b],
                                "nearest_branch_ratio": br,
                                "branch_ratio_error": err,
                            }
            row = {
                "k_low": lo,
                "k_high": hi,
                "ratio_high_low": ratio,
                "ratio_rational": rat,
                "ratio_rational_error": rat_err,
            }
            if branch_best:
                row.update(branch_best)
            rows.append(row)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows).sort_values(["branch_ratio_error", "ratio_rational_error"]).reset_index(drop=True)


def make_scan_dataset(event_df: pd.DataFrame, x_var: str, y_var: str, min_events: int) -> Optional[pd.DataFrame]:
    x_col = resolve_variable_column(event_df, x_var)
    y_col = resolve_variable_column(event_df, y_var)
    if x_col is None or y_col is None:
        return None
    x = pd.to_numeric(event_df[x_col], errors="coerce")
    y = pd.to_numeric(event_df[y_col], errors="coerce")
    y_base = y_col[:-len("_median")] if y_col.endswith("_median") else y_col
    w_col = f"{y_base}_width_16_84"
    if w_col in event_df.columns:
        width = pd.to_numeric(event_df[w_col], errors="coerce").to_numpy(float)
        w = 1.0 / np.maximum(width * width, EPS)
    else:
        w = np.ones(len(event_df), dtype=float)
    out = pd.DataFrame({"event": event_df["event"].astype(str), "x": x, "y": y, "w": w})
    out = out[np.isfinite(out["x"]) & np.isfinite(out["y"]) & np.isfinite(out["w"]) & (out["w"] > 0)].copy()
    if len(out) < min_events:
        return None
    if out["x"].std() <= 0 or out["y"].std() <= 0:
        return None
    return out.sort_values("x").reset_index(drop=True)


def choose_pairs(event_df: pd.DataFrame, requested: Optional[List[str]]) -> List[Tuple[str, str]]:
    if requested:
        pairs = []
        for item in requested:
            if ":" not in item:
                raise ValueError(f"--pair must be X:Y, got {item}")
            a, b = item.split(":", 1)
            pairs.append((a.strip(), b.strip()))
        return pairs

    x_candidates = []
    for v in SCAN_VARIABLES_POSITIVE:
        lv = f"log_{v}"
        if resolve_variable_column(event_df, lv) is not None:
            x_candidates.append(lv)
    y_candidates = [v for v in Y_VARIABLES if resolve_variable_column(event_df, v) is not None]
    return [(x, y) for x in x_candidates for y in y_candidates if x != y]


def build_event_table(files: List[str], out_dir: Path, posterior_choice: str, max_samples: Optional[int], use_weights: bool):
    rows = []
    audit_rows = []
    for i, f in enumerate(files, 1):
        print(f"[extract] {i}/{len(files)} {f}")
        try:
            row, tables = extract_event_summary(f, posterior_choice=posterior_choice, max_samples=max_samples, use_weights=use_weights)
            rows.append(row)
            audit_rows.extend(tables)
            if row.get("posterior_dataset"):
                print(f"  [posterior] {row['posterior_dataset']} n={row['posterior_n_samples']}")
            else:
                print("  [posterior] NONE")
        except Exception as exc:
            print(f"[warn] failed {f}: {exc}")
            rows.append({"event": event_name_from_file(f), "file": f, "extract_error": str(exc)})
    event_df = pd.DataFrame(rows)
    event_df = add_log_columns(event_df)
    event_df.to_csv(out_dir / "gwtc4_event_summary.csv", index=False)
    audit = pd.DataFrame(audit_rows)
    if not audit.empty:
        audit.to_csv(out_dir / "gwtc4_posterior_table_audit.csv", index=False)
    return event_df, audit


def run_scans(event_df: pd.DataFrame, out_dir: Path, pairs: List[Tuple[str, str]], k_grid: np.ndarray, degree: int, null_n: int, min_events: int, seed: int, skip_null: bool):
    scans_dir = out_dir / "scans"
    nulls_dir = out_dir / "nulls"
    scans_dir.mkdir(parents=True, exist_ok=True)
    nulls_dir.mkdir(parents=True, exist_ok=True)

    scan_rows, branch_rows, ratio_rows = [], [], []
    rng = np.random.default_rng(seed)

    for x_var, y_var in pairs:
        ds = make_scan_dataset(event_df, x_var, y_var, min_events=min_events)
        if ds is None:
            continue

        x = ds["x"].to_numpy(float)
        y = ds["y"].to_numpy(float)
        w = ds["w"].to_numpy(float)

        print("\n" + "=" * 80)
        print(f"[scan] {x_var} -> {y_var} N={len(ds)}")
        print("=" * 80)

        label = safe_name(f"{x_var}__{y_var}")
        scan_df, best, _, _ = harmonic_scan(x, y, w, k_grid, degree)
        scan_path = scans_dir / f"scan_{label}.csv"
        scan_df.to_csv(scan_path, index=False)

        ddom = delta_domain(x)
        best.update({
            "x_var": x_var,
            "y_var": y_var,
            "n_events": int(len(ds)),
            "x_min": float(np.min(x)),
            "x_max": float(np.max(x)),
            "delta_x": float(ddom),
            "n_obs": float(n_from_k(best["k_best"], ddom)),
            "scan_csv": str(scan_path),
        })

        if not skip_null and null_n > 0:
            null_seed = int(rng.integers(0, 2**31 - 1))
            null_max, p, count = null_scan(x, y, w, k_grid, degree, best["delta_chi2"], null_n, null_seed)
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

        bm = branch_grid_matches(best["k_best"], ddom).head(10)
        bm["x_var"] = x_var
        bm["y_var"] = y_var
        bm["scan_k_best"] = best["k_best"]
        bm["scan_delta_chi2"] = best["delta_chi2"]
        bm["scan_null_p"] = best["scan_null_p"]
        branch_rows.append(bm)

        rr = ratio_scan_results(scan_df)
        if not rr.empty:
            rr["x_var"] = x_var
            rr["y_var"] = y_var
            rr["scan_k_best"] = best["k_best"]
            rr["scan_delta_chi2"] = best["delta_chi2"]
            rr["scan_null_p"] = best["scan_null_p"]
            ratio_rows.append(rr)

        print(json.dumps({k: best[k] for k in ["x_var", "y_var", "n_events", "k_best", "delta_chi2", "n_obs", "scan_null_p", "tail_count_ge", "null_n"]}, indent=2))

    scans = pd.DataFrame(scan_rows)
    branches = pd.concat(branch_rows, ignore_index=True) if branch_rows else pd.DataFrame()
    ratios = pd.concat(ratio_rows, ignore_index=True) if ratio_rows else pd.DataFrame()

    if not scans.empty:
        scans = scans.sort_values(["scan_null_p", "delta_chi2"], ascending=[True, False], na_position="last").reset_index(drop=True)
        scans.to_csv(out_dir / "gwtc4_scan_results.csv", index=False)
    if not branches.empty:
        branches = branches.sort_values(["scan_null_p", "abs_k_error"], ascending=[True, True], na_position="last").reset_index(drop=True)
        branches.to_csv(out_dir / "gwtc4_branch_match_results.csv", index=False)
    if not ratios.empty:
        ratios = ratios.sort_values(["scan_null_p", "branch_ratio_error", "ratio_rational_error"], ascending=[True, True, True], na_position="last").reset_index(drop=True)
        ratios.to_csv(out_dir / "gwtc4_ratio_results.csv", index=False)

    return scans, branches, ratios


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", nargs="*", default=["*combined_PEDataRelease.hdf5"])
    ap.add_argument("--out-dir", default=OUTDIR_DEFAULT)
    ap.add_argument("--posterior-choice", default="preferred", choices=["preferred", "largest", "first"])
    ap.add_argument("--max-samples", type=int, default=None)
    ap.add_argument("--use-weights", action="store_true")

    ap.add_argument("--k-min", type=float, default=K_MIN_DEFAULT)
    ap.add_argument("--k-max", type=float, default=K_MAX_DEFAULT)
    ap.add_argument("--n-k", type=int, default=N_K_DEFAULT)
    ap.add_argument("--degree", type=int, default=BASELINE_DEGREE_DEFAULT)
    ap.add_argument("--null-n", type=int, default=NULL_N_DEFAULT)
    ap.add_argument("--skip-null", action="store_true")
    ap.add_argument("--seed", type=int, default=SEED_DEFAULT)
    ap.add_argument("--min-events", type=int, default=MIN_EVENTS_DEFAULT)
    ap.add_argument("--pair", action="append", default=None)
    ap.add_argument("--top-n", type=int, default=30)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = discover_files(args.input)
    files = [f for f in files if f.endswith(".hdf5") or f.endswith(".h5")]
    files = [f for f in files if "PESummaryTable" not in Path(f).name]
    if not files:
        raise FileNotFoundError(f"No event HDF5 files found for input {args.input}")

    print(f"[files] {len(files)} files")
    event_df, audit = build_event_table(files, out_dir, args.posterior_choice, args.max_samples, args.use_weights)

    usable_cols = [c for c in event_df.columns if c.endswith("_median")]
    print(f"[events] {len(event_df)} rows, median columns={len(usable_cols)}")
    if not usable_cols:
        raise RuntimeError("No posterior parameter columns extracted. Check gwtc4_posterior_table_audit.csv")

    pairs = choose_pairs(event_df, args.pair)
    print(f"[pairs] {len(pairs)} scan pairs")

    k_grid = np.linspace(args.k_min, args.k_max, args.n_k)
    scans, branches, ratios = run_scans(event_df, out_dir, pairs, k_grid, args.degree, args.null_n, args.min_events, args.seed, args.skip_null)

    summary = {
        "args": vars(args),
        "n_files": len(files),
        "n_events": int(len(event_df)),
        "n_posterior_tables_seen": int(len(audit)) if not audit.empty else 0,
        "n_scans": int(len(scans)),
        "n_branch_rows": int(len(branches)),
        "n_ratio_rows": int(len(ratios)),
        "delta_ell_lhc": DELTA_ELL_LHC,
        "top_scans": scans.head(args.top_n).to_dict(orient="records") if not scans.empty else [],
        "top_branch_matches": branches.head(args.top_n).to_dict(orient="records") if not branches.empty else [],
        "top_ratio_matches": ratios.head(args.top_n).to_dict(orient="records") if not ratios.empty else [],
    }
    with open(out_dir / "gwtc4_run_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n[done]")
    for fn in [
        "gwtc4_event_summary.csv",
        "gwtc4_posterior_table_audit.csv",
        "gwtc4_scan_results.csv",
        "gwtc4_branch_match_results.csv",
        "gwtc4_ratio_results.csv",
        "gwtc4_run_summary.json",
    ]:
        print("[save]", out_dir / fn)

    if not scans.empty:
        cols = ["x_var", "y_var", "n_events", "k_best", "delta_chi2", "n_obs", "scan_null_p", "tail_count_ge", "null_n"]
        print("\n[top scans]")
        print(scans.head(args.top_n)[cols].to_string(index=False))

    if not ratios.empty:
        cols = ["x_var", "y_var", "k_low", "k_high", "ratio_high_low", "ratio_rational", "nearest_branch_name", "nearest_branch_low", "nearest_branch_high", "branch_ratio_error", "scan_null_p"]
        print("\n[top ratios]")
        print(ratios.head(args.top_n)[cols].to_string(index=False))


if __name__ == "__main__":
    main()
