#!/usr/bin/env python3
"""Generate a dependence-preserving GWTC-4 population null matrix.

This script restores the missing outer/global-null layer for the historical
168-scan GWTC-4 subset analysis.

Historical manifest
-------------------
The recovered historical run contains 168 scans = 56 subset definitions x
3 final-spin pairs.  The original per-scan p-values were obtained from an
inner residual-permutation null independently for each scan.  Those p-values
are useful marginally, but their aggregate chi-square assumes more independence
than the overlapping selectors actually have.

Outer null used here
--------------------
One outer-null catalog is produced by permuting ``final_spin_median`` among the
events for which final spin is finite, while keeping event identity, masses,
chi_eff, q/eta, SNR, distance, missingness, and the event-specific final-spin
uncertainty/weight fixed.  The exact historical subset rules are then rebuilt.
This deliberately preserves:

- overlap among top-N / quantile selectors,
- reuse of the same events across the three x variables,
- selection on final spin itself (which is rebuilt after permutation),
- all non-final-spin catalog correlations and measurement precisions.

It tests the null that final-spin values are exchangeable across event labels
with respect to the historical mass/selector architecture.  It is therefore a
correlation-aware catalog-level permutation null, not a full LVK astrophysical
population model.

Calibration
-----------
For each historical scan definition, the observed scan-max Delta-chi2 and the
outer-null scan-max Delta-chi2 ensemble are converted to empirical upper-tail
p-values.  Null-row p-values use leave-one-out/exchangeable ranks.  This creates
``gwtc4_population_null_matrix.csv`` with one full 168-p-value vector per outer
null catalog and a separately calibrated observed 168-vector.

The resulting files can be passed directly to
``run_gwtc_population_global_null.py`` to calibrate the population-level
histogram statistic while preserving selector dependence.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from gwtc4_wct_subset_scan_compound import (
    EPS,
    FINAL_SPIN_PAIRS,
    add_derived_columns,
    build_subsets,
    make_poly_design,
    make_scan_dataset,
)


DYNAMIC_SUBSET_PATTERNS = (
    re.compile(r"^high_final_spin_"),
    re.compile(r"^low_final_spin_"),
    re.compile(r"^rank_top\d+_final_spin$"),
)


@dataclass
class PreparedScan:
    events: List[str]
    x: np.ndarray
    w: np.ndarray
    sw: np.ndarray
    c_resid: np.ndarray
    s_resid: np.ndarray
    aa: np.ndarray
    bb: np.ndarray
    cc: np.ndarray
    det: np.ndarray
    k_grid: np.ndarray


def is_dynamic_subset(name: str) -> bool:
    return any(p.search(str(name)) is not None for p in DYNAMIC_SUBSET_PATTERNS)


def validate_manifest(observed: pd.DataFrame, expected_n: int) -> pd.DataFrame:
    required = {"subset", "x_var", "y_var", "delta_chi2", "k_best", "n_events"}
    missing = required.difference(observed.columns)
    if missing:
        raise SystemExit(f"Observed manifest missing columns: {sorted(missing)}")
    if len(observed) != expected_n:
        raise SystemExit(f"Observed manifest has {len(observed)} rows; expected {expected_n}")

    out = observed.copy().reset_index(drop=True)
    out["scan_id"] = [f"scan_{i:03d}" for i in range(len(out))]

    pairs = set(zip(out["x_var"].astype(str), out["y_var"].astype(str)))
    expected_pairs = set(FINAL_SPIN_PAIRS)
    if pairs != expected_pairs:
        raise SystemExit(
            "Historical manifest pair set does not match FINAL_SPIN_PAIRS. "
            f"observed={sorted(pairs)!r}, expected={sorted(expected_pairs)!r}"
        )

    if out[["subset", "x_var", "y_var"]].duplicated().any():
        dup = out.loc[
            out[["subset", "x_var", "y_var"]].duplicated(keep=False),
            ["subset", "x_var", "y_var"],
        ]
        raise SystemExit(f"Duplicate scan identities in manifest:\n{dup.to_string(index=False)}")
    return out


def prepare_scan(
    subset_df: pd.DataFrame,
    x_var: str,
    y_var: str,
    min_events: int,
    k_grid: np.ndarray,
    degree: int,
) -> PreparedScan:
    ds = make_scan_dataset(subset_df, x_var, y_var, min_events)
    if ds is None:
        raise ValueError(f"Could not build scan dataset for {x_var}->{y_var}")

    x = ds["x"].to_numpy(float)
    w = ds["w"].to_numpy(float)
    sw = np.sqrt(np.maximum(w, EPS))

    x0w = make_poly_design(x, degree) * sw[:, None]
    c = np.cos(np.outer(x, k_grid)) * sw[:, None]
    s = np.sin(np.outer(x, k_grid)) * sw[:, None]

    # Frisch-Waugh-Lovell: residualize the two harmonic columns against the
    # polynomial baseline in weighted space.  The resulting Delta-chi2 is
    # algebraically identical to refitting [X0, cos(kx), sin(kx)] at every k.
    c_resid = c - x0w @ np.linalg.lstsq(x0w, c, rcond=None)[0]
    s_resid = s - x0w @ np.linalg.lstsq(x0w, s, rcond=None)[0]

    aa = np.sum(c_resid * c_resid, axis=0)
    bb = np.sum(c_resid * s_resid, axis=0)
    cc = np.sum(s_resid * s_resid, axis=0)
    det = aa * cc - bb * bb

    return PreparedScan(
        events=ds["event"].astype(str).tolist(),
        x=x,
        w=w,
        sw=sw,
        c_resid=c_resid,
        s_resid=s_resid,
        aa=aa,
        bb=bb,
        cc=cc,
        det=det,
        k_grid=k_grid,
    )


def score_prepared_one(prep: PreparedScan, y: np.ndarray) -> Tuple[float, float]:
    y = np.asarray(y, dtype=float)
    if y.shape != prep.sw.shape:
        raise ValueError(f"y shape {y.shape} != expected {prep.sw.shape}")

    yw = prep.sw * y
    u = prep.c_resid.T @ yw
    v = prep.s_resid.T @ yw
    numer = prep.cc * u * u - 2.0 * prep.bb * u * v + prep.aa * v * v
    good = prep.det > 1e-14
    delta = np.zeros_like(prep.det)
    delta[good] = numer[good] / prep.det[good]
    delta = np.maximum(delta, 0.0)
    idx = int(np.argmax(delta))
    return float(prep.k_grid[idx]), float(delta[idx])


def score_prepared_many(prep: PreparedScan, y_matrix: np.ndarray, batch_size: int = 256) -> np.ndarray:
    """Score many y-vectors sharing the same x, weights and event membership."""
    y_matrix = np.asarray(y_matrix, dtype=float)
    if y_matrix.ndim != 2 or y_matrix.shape[1] != len(prep.events):
        raise ValueError(
            f"Expected y_matrix shape (M,{len(prep.events)}), got {y_matrix.shape}"
        )

    out = np.empty(y_matrix.shape[0], dtype=float)
    good = prep.det > 1e-14

    for start in range(0, y_matrix.shape[0], batch_size):
        stop = min(start + batch_size, y_matrix.shape[0])
        yw = prep.sw[:, None] * y_matrix[start:stop].T
        u = prep.c_resid.T @ yw
        v = prep.s_resid.T @ yw
        numer = (
            prep.cc[:, None] * u * u
            - 2.0 * prep.bb[:, None] * u * v
            + prep.aa[:, None] * v * v
        )
        delta = np.zeros_like(numer)
        delta[good, :] = numer[good, :] / prep.det[good, None]
        delta = np.maximum(delta, 0.0)
        out[start:stop] = np.max(delta, axis=0)
    return out


def map_subsets(
    event_df: pd.DataFrame,
    min_events: int,
    top_ns: List[int],
) -> Dict[str, pd.DataFrame]:
    subsets = build_subsets(event_df, min_events=min_events, top_ns=top_ns, custom_cuts=None)
    return {str(name): sub.copy() for name, sub in subsets}


def event_value_map(df: pd.DataFrame, column: str) -> Dict[str, float]:
    if column not in df.columns:
        raise SystemExit(f"Required event-summary column not found: {column}")
    return dict(zip(df["event"].astype(str), pd.to_numeric(df[column], errors="coerce")))


def y_for_events(values: Dict[str, float], events: Iterable[str]) -> np.ndarray:
    arr = np.array([values.get(str(e), np.nan) for e in events], dtype=float)
    if not np.all(np.isfinite(arr)):
        bad = [str(e) for e, v in zip(events, arr) if not np.isfinite(v)]
        raise ValueError(f"Non-finite final spin for prepared events: {bad[:5]}")
    return arr


def make_outer_permutations(event_df: pd.DataFrame, outer_n: int, seed: int) -> Tuple[np.ndarray, np.ndarray]:
    y = pd.to_numeric(event_df["final_spin_median"], errors="coerce").to_numpy(float)
    finite_idx = np.flatnonzero(np.isfinite(y))
    if finite_idx.size < 2:
        raise SystemExit("Need at least two finite final_spin_median values")

    rng = np.random.default_rng(seed)
    y_perm = np.tile(y, (outer_n, 1))
    finite_values = y[finite_idx].copy()
    for i in range(outer_n):
        y_perm[i, finite_idx] = finite_values[rng.permutation(finite_idx.size)]
    return y_perm, finite_idx


def upper_rank_p_observed(null_delta: np.ndarray, observed_delta: np.ndarray) -> np.ndarray:
    return (1.0 + np.sum(null_delta >= observed_delta[None, :], axis=0)) / (1.0 + null_delta.shape[0])


def leave_one_out_rank_p(null_delta: np.ndarray) -> np.ndarray:
    """Exchangeable upper-tail ranks for each null row, column by column."""
    m, n = null_delta.shape
    p = np.empty((m, n), dtype=float)
    for j in range(n):
        # rankdata(-x, method='max') equals count(values >= x), including self.
        p[:, j] = rankdata(-null_delta[:, j], method="max") / float(m)
    return p


def main() -> None:
    ap = argparse.ArgumentParser(description="Generate correlated GWTC-4 168-scan outer-null matrix")
    ap.add_argument(
        "--event-summary",
        default="outputs_gwtc4_wct_subsets/gwtc4_subset_event_summary.csv",
        help="Recovered historical event-summary CSV",
    )
    ap.add_argument(
        "--observed",
        default="tables/gwtc4_population_observed.csv",
        help="Recovered historical 168-row scan manifest/results",
    )
    ap.add_argument("--expected-n", type=int, default=168)
    ap.add_argument("--outer-n", type=int, default=200)
    ap.add_argument("--seed", type=int, default=20260831)
    ap.add_argument("--k-min", type=float, default=0.5)
    ap.add_argument("--k-max", type=float, default=80.0)
    ap.add_argument("--n-k", type=int, default=4000)
    ap.add_argument("--degree", type=int, default=2)
    ap.add_argument("--min-events", type=int, default=12)
    ap.add_argument("--top-ns", default="12,16,20,24,32")
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument("--verify-only", action="store_true")
    ap.add_argument(
        "--output-matrix", default="tables/gwtc4_population_null_matrix.csv"
    )
    ap.add_argument(
        "--output-observed", default="tables/gwtc4_population_observed_outercal.csv"
    )
    ap.add_argument(
        "--output-delta", default="tables/gwtc4_population_null_delta_matrix.csv"
    )
    ap.add_argument(
        "--metadata", default="tables/gwtc4_population_null_metadata.json"
    )
    args = ap.parse_args()

    if args.outer_n < 2 and not args.verify_only:
        raise SystemExit("--outer-n must be at least 2")
    if args.n_k < 2:
        raise SystemExit("--n-k must be at least 2")

    event_path = Path(args.event_summary)
    observed_path = Path(args.observed)
    if not event_path.exists():
        raise SystemExit(f"Historical event-summary file not found: {event_path}")
    if not observed_path.exists():
        raise SystemExit(f"Historical observed manifest not found: {observed_path}")

    event_df = add_derived_columns(pd.read_csv(event_path))
    if "event" not in event_df.columns:
        raise SystemExit("Event summary lacks 'event' column")
    if event_df["event"].astype(str).duplicated().any():
        raise SystemExit("Event summary contains duplicate event names")

    observed = validate_manifest(pd.read_csv(observed_path), args.expected_n)
    top_ns = [int(x.strip()) for x in str(args.top_ns).split(",") if x.strip()]
    k_grid = np.linspace(args.k_min, args.k_max, args.n_k)

    subset_map = map_subsets(event_df, args.min_events, top_ns)
    manifest_subsets = set(observed["subset"].astype(str))
    missing_subsets = sorted(manifest_subsets.difference(subset_map))
    if missing_subsets:
        raise SystemExit(
            "Recovered event summary does not reproduce historical subset names. "
            f"Missing: {missing_subsets[:10]}"
        )

    print("Historical manifest verification")
    print("--------------------------------")
    print(f"Events                  : {len(event_df)}")
    print(f"Observed scans          : {len(observed)}")
    print(f"Unique subsets          : {len(manifest_subsets)}")
    print(f"Pair definitions        : {len(set(zip(observed.x_var, observed.y_var)))}")
    print(f"Dynamic final-spin sets : {sum(is_dynamic_subset(s) for s in manifest_subsets)}")

    # Recompute the observed scan-max statistics from the recovered event summary
    # before generating any nulls.  This guards against silently using a different
    # grid, selector set, event table, weighting rule, or baseline degree.
    observed_recalc_k = np.empty(len(observed), dtype=float)
    observed_recalc_delta = np.empty(len(observed), dtype=float)
    prepared_fixed: Dict[int, PreparedScan] = {}

    original_y_map = event_value_map(event_df, "final_spin_median")

    for i, row in observed.iterrows():
        subset_name = str(row["subset"])
        prep = prepare_scan(
            subset_map[subset_name],
            str(row["x_var"]),
            str(row["y_var"]),
            args.min_events,
            k_grid,
            args.degree,
        )
        y = y_for_events(original_y_map, prep.events)
        k_best, delta = score_prepared_one(prep, y)
        observed_recalc_k[i] = k_best
        observed_recalc_delta[i] = delta
        if not is_dynamic_subset(subset_name):
            prepared_fixed[i] = prep

    k_err = np.abs(observed_recalc_k - pd.to_numeric(observed["k_best"]).to_numpy(float))
    d_err = np.abs(observed_recalc_delta - pd.to_numeric(observed["delta_chi2"]).to_numpy(float))
    max_k_err = float(np.max(k_err))
    max_d_err = float(np.max(d_err))

    print(f"Max |k_recalc-k_saved|   : {max_k_err:.6g}")
    print(f"Max |D_recalc-D_saved|   : {max_d_err:.6g}")

    # Exact recovered settings should reproduce to floating precision.  A loose
    # 1e-7 Delta tolerance allows platform-level BLAS roundoff but catches a
    # changed grid/weight/selector definition immediately.
    if max_k_err > 1e-10 or max_d_err > 1e-7:
        raise SystemExit(
            "Observed verification failed: recovered event summary/settings do not "
            "reproduce the historical 168 scan statistics."
        )

    print("Observed reproduction     : PASS")

    metadata = {
        "schema": "GWTC4_POPULATION_OUTER_NULL_V1",
        "null_type": "final_spin_event_label_permutation",
        "interpretation": (
            "Permute finite final_spin_median values across event labels while "
            "holding event-specific uncertainty weights and all non-final-spin "
            "catalog quantities fixed; rebuild all historical selectors."
        ),
        "event_summary": str(event_path),
        "observed_manifest": str(observed_path),
        "n_events": int(len(event_df)),
        "n_scans": int(len(observed)),
        "n_unique_subsets": int(len(manifest_subsets)),
        "n_dynamic_final_spin_subsets": int(sum(is_dynamic_subset(s) for s in manifest_subsets)),
        "pairs": sorted([list(x) for x in set(zip(observed.x_var.astype(str), observed.y_var.astype(str)))]),
        "k_min": float(args.k_min),
        "k_max": float(args.k_max),
        "n_k": int(args.n_k),
        "degree": int(args.degree),
        "min_events": int(args.min_events),
        "top_ns": top_ns,
        "seed": int(args.seed),
        "outer_n": 0 if args.verify_only else int(args.outer_n),
        "observed_reproduction_max_abs_k_error": max_k_err,
        "observed_reproduction_max_abs_delta_error": max_d_err,
        "observed_reproduction_pass": True,
    }

    meta_path = Path(args.metadata)
    meta_path.parent.mkdir(parents=True, exist_ok=True)

    if args.verify_only:
        meta_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
        print(f"Wrote metadata            : {meta_path}")
        return

    # Create the shared outer-null catalogs once at the event level.
    y_perm, finite_idx = make_outer_permutations(event_df, args.outer_n, args.seed)
    metadata["n_finite_final_spin"] = int(len(finite_idx))

    event_names = event_df["event"].astype(str).tolist()
    event_index = {e: i for i, e in enumerate(event_names)}
    null_delta = np.empty((args.outer_n, len(observed)), dtype=float)

    fixed_indices = [i for i, row in observed.iterrows() if not is_dynamic_subset(str(row["subset"]))]
    dynamic_indices = [i for i, row in observed.iterrows() if is_dynamic_subset(str(row["subset"]))]

    print("\nOuter-null scoring")
    print("------------------")
    print(f"Outer catalogs           : {args.outer_n}")
    print(f"Fixed-membership scans   : {len(fixed_indices)}")
    print(f"Dynamic-selector scans   : {len(dynamic_indices)}")

    # Fixed-membership scans are scored for all outer catalogs in BLAS batches.
    for pos, i in enumerate(fixed_indices, 1):
        prep = prepared_fixed[i]
        cols = [event_index[e] for e in prep.events]
        y_many = y_perm[:, cols]
        null_delta[:, i] = score_prepared_many(prep, y_many, batch_size=args.batch_size)
        if pos % 20 == 0 or pos == len(fixed_indices):
            print(f"[fixed] {pos}/{len(fixed_indices)}")

    # Final-spin-selected subsets change membership after each permutation, so
    # rebuild those subset definitions for every outer catalog.
    base_y = pd.to_numeric(event_df["final_spin_median"], errors="coerce").to_numpy(float)
    for outer_i in range(args.outer_n):
        null_df = event_df.copy()
        null_df["final_spin_median"] = y_perm[outer_i]
        null_df = add_derived_columns(null_df)
        dynamic_map = map_subsets(null_df, args.min_events, top_ns)
        y_map = event_value_map(null_df, "final_spin_median")

        for i in dynamic_indices:
            row = observed.iloc[i]
            subset_name = str(row["subset"])
            if subset_name not in dynamic_map:
                raise SystemExit(
                    f"Outer catalog {outer_i}: dynamic subset disappeared: {subset_name}"
                )
            prep = prepare_scan(
                dynamic_map[subset_name],
                str(row["x_var"]),
                str(row["y_var"]),
                args.min_events,
                k_grid,
                args.degree,
            )
            y = y_for_events(y_map, prep.events)
            _, delta = score_prepared_one(prep, y)
            null_delta[outer_i, i] = delta

        if (outer_i + 1) % 25 == 0 or outer_i + 1 == args.outer_n:
            print(f"[dynamic] {outer_i + 1}/{args.outer_n}")

    if not np.all(np.isfinite(null_delta)):
        raise SystemExit("Generated null Delta matrix contains non-finite values")

    observed_delta = pd.to_numeric(observed["delta_chi2"], errors="coerce").to_numpy(float)
    observed_outer_p = upper_rank_p_observed(null_delta, observed_delta)
    null_p = leave_one_out_rank_p(null_delta)

    scan_columns = observed["scan_id"].astype(str).tolist()
    matrix_df = pd.DataFrame(null_p, columns=scan_columns)
    delta_df = pd.DataFrame(null_delta, columns=scan_columns)

    observed_outer = observed[
        ["scan_id", "subset", "x_var", "y_var", "n_events", "k_best", "delta_chi2"]
    ].copy()
    observed_outer["outer_global_p"] = observed_outer_p

    matrix_path = Path(args.output_matrix)
    observed_out_path = Path(args.output_observed)
    delta_path = Path(args.output_delta)
    for p in [matrix_path, observed_out_path, delta_path]:
        p.parent.mkdir(parents=True, exist_ok=True)

    matrix_df.to_csv(matrix_path, index=False)
    delta_df.to_csv(delta_path, index=False)
    observed_outer.to_csv(observed_out_path, index=False)

    metadata.update(
        {
            "output_null_p_matrix": str(matrix_path),
            "output_null_delta_matrix": str(delta_path),
            "output_observed_outercal": str(observed_out_path),
            "observed_outer_p_min": float(np.min(observed_outer_p)),
            "observed_outer_p_median": float(np.median(observed_outer_p)),
            "observed_outer_p_lt_0p05": int(np.sum(observed_outer_p < 0.05)),
            "observed_outer_p_lt_0p10": int(np.sum(observed_outer_p < 0.10)),
            "outer_per_scan_resolution_floor": float(1.0 / (1.0 + args.outer_n)),
            "null_row_rank_resolution_floor": float(1.0 / args.outer_n),
        }
    )
    meta_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    print("\nGenerated")
    print("---------")
    print(f"Observed outer-cal vector : {observed_out_path}")
    print(f"Null p-value matrix       : {matrix_path}")
    print(f"Null Delta matrix         : {delta_path}")
    print(f"Metadata                  : {meta_path}")
    print(f"Per-scan p floor          : {1.0 / (1.0 + args.outer_n):.8g}")
    print(f"Observed p<0.05 count     : {int(np.sum(observed_outer_p < 0.05))}")
    print(f"Observed p<0.10 count     : {int(np.sum(observed_outer_p < 0.10))}")
    print("\nNext command:")
    print(
        "python scripts/run_gwtc_population_global_null.py "
        f"--observed {observed_out_path} --p-column outer_global_p "
        f"--null-matrix {matrix_path} --output tables/gwtc4_population_global_result.json"
    )


if __name__ == "__main__":
    main()
