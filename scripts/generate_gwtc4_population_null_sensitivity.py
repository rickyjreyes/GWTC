#!/usr/bin/env python3
"""GWTC-4 population-null sensitivity study.

Keep the recovered historical 168-scan statistic and selector architecture
fixed while changing only the catalog-level null model. This is intentionally
separate from ``generate_gwtc4_population_null_matrix.py`` so the established
full event-label permutation result remains frozen/reproducible.

Null models
-----------
label
    Full final-spin event-label permutation. Included as a cross-check.

mass_stratified
    Permute final-spin medians only within quantile strata of a declared mass
    coordinate. This preserves broad mass association while destroying finer
    event-level ordering.

mass_precision_stratified
    Permute final-spin medians within joint mass-quantile and measurement-
    precision-quantile cells. Event-specific historical weights remain fixed.

mass_residual
    Fit a smooth weighted polynomial relation between a declared mass
    coordinate and atanh(final_spin), permute fitted residuals across events,
    then transform back with tanh. This preserves a smooth broad mass-spin
    relation while destroying fine residual ordering and keeps generated spins
    inside (-1, 1).
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

from generate_gwtc4_population_null_matrix import (
    PreparedScan,
    event_value_map,
    is_dynamic_subset,
    leave_one_out_rank_p,
    map_subsets,
    prepare_scan,
    score_prepared_many,
    score_prepared_one,
    upper_rank_p_observed,
    validate_manifest,
    y_for_events,
)
from gwtc4_wct_subset_scan_compound import (
    EPS,
    add_derived_columns,
    make_poly_design,
    resolve_variable_column,
)


NULL_MODELS = (
    "label",
    "mass_stratified",
    "mass_precision_stratified",
    "mass_residual",
)


def numeric_column(df: pd.DataFrame, variable: str) -> Tuple[str, np.ndarray]:
    col = resolve_variable_column(df, variable)
    if col is None and variable in df.columns:
        col = variable
    if col is None:
        raise SystemExit(f"Could not resolve variable/column: {variable}")
    return col, pd.to_numeric(df[col], errors="coerce").to_numpy(float)


def finite_spin(event_df: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
    y = pd.to_numeric(event_df["final_spin_median"], errors="coerce").to_numpy(float)
    idx = np.flatnonzero(np.isfinite(y))
    if idx.size < 2:
        raise SystemExit("Need at least two finite final_spin_median values")
    return y, idx


def quantile_codes(
    values: np.ndarray, indices: np.ndarray, bins: int
) -> Tuple[np.ndarray, Dict[str, object]]:
    if bins < 2:
        raise SystemExit("Quantile stratification requires at least 2 bins")
    vals = np.asarray(values, dtype=float)[indices]
    if not np.all(np.isfinite(vals)):
        raise SystemExit(
            "Stratification variable is non-finite for one or more finite-spin events"
        )
    unique = np.unique(vals)
    if unique.size < 2:
        raise SystemExit("Stratification variable has fewer than two unique values")

    q = min(int(bins), int(unique.size))
    codes = pd.qcut(
        pd.Series(vals), q=q, labels=False, duplicates="drop"
    ).to_numpy()
    if np.any(pd.isna(codes)):
        raise SystemExit("Could not assign all finite-spin events to quantile strata")
    codes = codes.astype(int)
    counts = {str(int(k)): int(np.sum(codes == k)) for k in np.unique(codes)}
    return codes, {
        "requested_bins": int(bins),
        "realized_bins": int(len(counts)),
        "cell_sizes": counts,
    }


def permute_within_groups(
    y: np.ndarray,
    finite_idx: np.ndarray,
    groups: np.ndarray,
    outer_n: int,
    rng: np.random.Generator,
) -> np.ndarray:
    y_perm = np.tile(y, (outer_n, 1))
    groups = np.asarray(groups)
    if groups.shape != finite_idx.shape:
        raise ValueError("groups and finite_idx must have the same shape")

    for outer_i in range(outer_n):
        for group in np.unique(groups):
            local = finite_idx[groups == group]
            if local.size > 1:
                y_perm[outer_i, local] = y[local][rng.permutation(local.size)]
    return y_perm


def make_label_nulls(
    event_df: pd.DataFrame, outer_n: int, seed: int
) -> Tuple[np.ndarray, np.ndarray, Dict[str, object]]:
    y, finite_idx = finite_spin(event_df)
    rng = np.random.default_rng(seed)
    groups = np.zeros(finite_idx.size, dtype=int)
    y_perm = permute_within_groups(y, finite_idx, groups, outer_n, rng)
    return y_perm, finite_idx, {
        "null_type": "final_spin_event_label_permutation",
        "interpretation": (
            "Permute finite final_spin_median values across all finite-spin event "
            "labels; event-specific weights and non-final-spin quantities remain fixed."
        ),
    }


def make_mass_stratified_nulls(
    event_df: pd.DataFrame,
    outer_n: int,
    seed: int,
    mass_var: str,
    mass_bins: int,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, object]]:
    y, finite_idx = finite_spin(event_df)
    mass_col, mass = numeric_column(event_df, mass_var)
    mass_codes, mass_meta = quantile_codes(mass, finite_idx, mass_bins)
    rng = np.random.default_rng(seed)
    y_perm = permute_within_groups(y, finite_idx, mass_codes, outer_n, rng)
    return y_perm, finite_idx, {
        "null_type": "final_spin_mass_stratified_permutation",
        "interpretation": (
            "Permute final_spin_median only within quantile strata of a fixed mass "
            "coordinate, preserving broad mass-spin association while destroying "
            "event-level ordering."
        ),
        "mass_variable": mass_var,
        "mass_column": mass_col,
        "mass_strata": mass_meta,
    }


def make_mass_precision_stratified_nulls(
    event_df: pd.DataFrame,
    outer_n: int,
    seed: int,
    mass_var: str,
    mass_bins: int,
    precision_bins: int,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, object]]:
    y, finite_idx = finite_spin(event_df)
    mass_col, mass = numeric_column(event_df, mass_var)
    width_col, width = numeric_column(event_df, "final_spin_width_16_84")
    if np.any(~np.isfinite(width[finite_idx]) | (width[finite_idx] <= 0)):
        raise SystemExit(
            "mass_precision_stratified requires finite positive final-spin widths"
        )

    mass_codes, mass_meta = quantile_codes(mass, finite_idx, mass_bins)
    width_codes, width_meta = quantile_codes(width, finite_idx, precision_bins)
    n_width = int(np.max(width_codes)) + 1
    joint = mass_codes * n_width + width_codes
    joint_sizes = {str(int(k)): int(np.sum(joint == k)) for k in np.unique(joint)}

    rng = np.random.default_rng(seed)
    y_perm = permute_within_groups(y, finite_idx, joint, outer_n, rng)
    return y_perm, finite_idx, {
        "null_type": "final_spin_mass_precision_stratified_permutation",
        "interpretation": (
            "Permute final_spin_median only within joint mass-quantile and final-spin-"
            "uncertainty-quantile cells. Event-specific uncertainty weights remain "
            "fixed, so donor and recipient events have broadly similar mass and precision."
        ),
        "mass_variable": mass_var,
        "mass_column": mass_col,
        "precision_column": width_col,
        "mass_strata": mass_meta,
        "precision_strata": width_meta,
        "joint_cell_sizes": joint_sizes,
        "singleton_joint_cells": int(sum(v == 1 for v in joint_sizes.values())),
    }


def make_mass_residual_nulls(
    event_df: pd.DataFrame,
    outer_n: int,
    seed: int,
    mass_var: str,
    residual_degree: int,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, object]]:
    if residual_degree < 0:
        raise SystemExit("--residual-degree must be >= 0")

    y, finite_idx = finite_spin(event_df)
    mass_col, mass = numeric_column(event_df, mass_var)
    if not np.all(np.isfinite(mass[finite_idx])):
        raise SystemExit(
            "mass_residual requires finite mass coordinate for all finite-spin events"
        )
    if np.any(np.abs(y[finite_idx]) >= 1.0):
        raise SystemExit("mass_residual requires |final_spin_median| < 1")

    width_col, width = numeric_column(event_df, "final_spin_width_16_84")
    if np.any(~np.isfinite(width[finite_idx]) | (width[finite_idx] <= 0)):
        raise SystemExit("mass_residual requires finite positive final-spin widths")

    x = mass[finite_idx]
    z = np.arctanh(y[finite_idx])
    w = 1.0 / np.maximum(width[finite_idx] ** 2, EPS)
    X = make_poly_design(x, residual_degree)
    sw = np.sqrt(np.maximum(w, EPS))
    beta, *_ = np.linalg.lstsq(X * sw[:, None], z * sw, rcond=None)
    fitted = X @ beta
    resid = z - fitted

    rng = np.random.default_rng(seed)
    y_perm = np.tile(y, (outer_n, 1))
    for outer_i in range(outer_n):
        z_null = fitted + resid[rng.permutation(resid.size)]
        y_perm[outer_i, finite_idx] = np.tanh(z_null)

    return y_perm, finite_idx, {
        "null_type": "smooth_mass_spin_residual_permutation",
        "interpretation": (
            "Fit a weighted polynomial baseline in atanh(final_spin) versus the "
            "declared mass coordinate, permute residuals across events, and transform "
            "back with tanh. This preserves a smooth broad mass-spin relation while "
            "destroying fine residual ordering."
        ),
        "mass_variable": mass_var,
        "mass_column": mass_col,
        "precision_column": width_col,
        "residual_degree": int(residual_degree),
        "baseline_coefficients": [float(v) for v in beta],
        "residual_mean": float(np.mean(resid)),
        "residual_std": float(np.std(resid)),
        "generated_final_spin_min": float(np.nanmin(y_perm[:, finite_idx])),
        "generated_final_spin_max": float(np.nanmax(y_perm[:, finite_idx])),
    }


def generate_outer_y(
    event_df: pd.DataFrame,
    outer_n: int,
    seed: int,
    null_model: str,
    mass_var: str,
    mass_bins: int,
    precision_bins: int,
    residual_degree: int,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, object]]:
    if null_model == "label":
        return make_label_nulls(event_df, outer_n, seed)
    if null_model == "mass_stratified":
        return make_mass_stratified_nulls(
            event_df, outer_n, seed, mass_var, mass_bins
        )
    if null_model == "mass_precision_stratified":
        return make_mass_precision_stratified_nulls(
            event_df, outer_n, seed, mass_var, mass_bins, precision_bins
        )
    if null_model == "mass_residual":
        return make_mass_residual_nulls(
            event_df, outer_n, seed, mass_var, residual_degree
        )
    raise ValueError(f"Unknown null model: {null_model}")


def resolve_dynamic_subset(
    subset_map: Dict[str, pd.DataFrame], historical_name: str
) -> pd.DataFrame:
    """Resolve dynamic quantile subsets when null thresholds change their names."""
    if historical_name in subset_map:
        return subset_map[historical_name]

    prefixes = []
    if historical_name.startswith("high_final_spin_final_spin_ge_q0.75_"):
        prefixes.append("high_final_spin_final_spin_ge_q0.75_")
    if historical_name.startswith("low_final_spin_final_spin_le_q0.25_"):
        prefixes.append("low_final_spin_final_spin_le_q0.25_")

    for prefix in prefixes:
        matches = [name for name in subset_map if name.startswith(prefix)]
        if len(matches) == 1:
            return subset_map[matches[0]]
        if len(matches) > 1:
            raise SystemExit(
                f"Ambiguous dynamic subset resolution for {historical_name}: {matches}"
            )

    raise SystemExit(f"Dynamic subset not found under null catalog: {historical_name}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Generate alternative correlation-aware GWTC-4 population nulls"
    )
    ap.add_argument(
        "--event-summary",
        default="outputs_gwtc4_wct_subsets/gwtc4_subset_event_summary.csv",
    )
    ap.add_argument("--observed", default="tables/gwtc4_population_observed.csv")
    ap.add_argument("--expected-n", type=int, default=168)
    ap.add_argument("--outer-n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260902)
    ap.add_argument("--null-model", choices=NULL_MODELS, required=True)
    ap.add_argument("--mass-var", default="log_total_mass_source")
    ap.add_argument("--mass-bins", type=int, default=4)
    ap.add_argument("--precision-bins", type=int, default=2)
    ap.add_argument("--residual-degree", type=int, default=2)
    ap.add_argument("--k-min", type=float, default=0.5)
    ap.add_argument("--k-max", type=float, default=80.0)
    ap.add_argument("--n-k", type=int, default=4000)
    ap.add_argument("--degree", type=int, default=2)
    ap.add_argument("--min-events", type=int, default=12)
    ap.add_argument("--top-ns", default="12,16,20,24,32")
    ap.add_argument("--batch-size", type=int, default=128)
    ap.add_argument(
        "--output-prefix",
        default=None,
        help="Prefix without suffix; defaults to tables/gwtc4_population_<model>_<N>",
    )
    args = ap.parse_args()

    if args.outer_n < 2:
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
    print(
        f"Dynamic final-spin sets : {sum(is_dynamic_subset(s) for s in manifest_subsets)}"
    )

    original_y_map = event_value_map(event_df, "final_spin_median")
    observed_recalc_k = np.empty(len(observed), dtype=float)
    observed_recalc_delta = np.empty(len(observed), dtype=float)
    prepared_fixed: Dict[int, PreparedScan] = {}

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
        y_obs = y_for_events(original_y_map, prep.events)
        k_best, delta = score_prepared_one(prep, y_obs)
        observed_recalc_k[i] = k_best
        observed_recalc_delta[i] = delta
        if not is_dynamic_subset(subset_name):
            prepared_fixed[i] = prep

    k_saved = pd.to_numeric(observed["k_best"], errors="coerce").to_numpy(float)
    d_saved = pd.to_numeric(observed["delta_chi2"], errors="coerce").to_numpy(float)
    max_k_err = float(np.max(np.abs(observed_recalc_k - k_saved)))
    max_d_err = float(np.max(np.abs(observed_recalc_delta - d_saved)))
    print(f"Max |k_recalc-k_saved|   : {max_k_err:.6g}")
    print(f"Max |D_recalc-D_saved|   : {max_d_err:.6g}")
    if max_k_err > 1e-10 or max_d_err > 1e-7:
        raise SystemExit(
            "Observed historical reproduction failed; aborting null sensitivity run"
        )
    print("Observed reproduction     : PASS")

    y_perm, finite_idx, null_meta = generate_outer_y(
        event_df=event_df,
        outer_n=args.outer_n,
        seed=args.seed,
        null_model=args.null_model,
        mass_var=args.mass_var,
        mass_bins=args.mass_bins,
        precision_bins=args.precision_bins,
        residual_degree=args.residual_degree,
    )
    if y_perm.shape != (args.outer_n, len(event_df)):
        raise SystemExit(f"Unexpected outer-y shape: {y_perm.shape}")
    if not np.all(np.isfinite(y_perm[:, finite_idx])):
        raise SystemExit("Generated null catalogs contain non-finite final spins")

    event_names = event_df["event"].astype(str).tolist()
    event_index = {e: i for i, e in enumerate(event_names)}
    null_delta = np.empty((args.outer_n, len(observed)), dtype=float)
    fixed_indices = [
        i
        for i, row in observed.iterrows()
        if not is_dynamic_subset(str(row["subset"]))
    ]
    dynamic_indices = [
        i for i, row in observed.iterrows() if is_dynamic_subset(str(row["subset"]))
    ]

    print("\nNull sensitivity scoring")
    print("------------------------")
    print(f"Null model               : {args.null_model}")
    print(f"Outer catalogs           : {args.outer_n}")
    print(f"Fixed-membership scans   : {len(fixed_indices)}")
    print(f"Dynamic-selector scans   : {len(dynamic_indices)}")

    for pos, i in enumerate(fixed_indices, 1):
        prep = prepared_fixed[i]
        cols = [event_index[e] for e in prep.events]
        null_delta[:, i] = score_prepared_many(
            prep, y_perm[:, cols], batch_size=args.batch_size
        )
        if pos % 20 == 0 or pos == len(fixed_indices):
            print(f"[fixed] {pos}/{len(fixed_indices)}")

    for outer_i in range(args.outer_n):
        null_df = event_df.copy()
        null_df["final_spin_median"] = y_perm[outer_i]
        null_df = add_derived_columns(null_df)
        dynamic_map = map_subsets(null_df, args.min_events, top_ns)
        y_map = event_value_map(null_df, "final_spin_median")

        for i in dynamic_indices:
            row = observed.iloc[i]
            subset_name = str(row["subset"])
            dynamic_subset = resolve_dynamic_subset(dynamic_map, subset_name)
            prep = prepare_scan(
                dynamic_subset,
                str(row["x_var"]),
                str(row["y_var"]),
                args.min_events,
                k_grid,
                args.degree,
            )
            y_null = y_for_events(y_map, prep.events)
            _, delta = score_prepared_one(prep, y_null)
            null_delta[outer_i, i] = delta

        if (outer_i + 1) % 25 == 0 or outer_i + 1 == args.outer_n:
            print(f"[dynamic] {outer_i + 1}/{args.outer_n}")

    if not np.all(np.isfinite(null_delta)):
        raise SystemExit("Generated null Delta matrix contains non-finite values")

    observed_outer_p = upper_rank_p_observed(null_delta, d_saved)
    null_p = leave_one_out_rank_p(null_delta)

    scan_columns = observed["scan_id"].astype(str).tolist()
    matrix_df = pd.DataFrame(null_p, columns=scan_columns)
    delta_df = pd.DataFrame(null_delta, columns=scan_columns)
    observed_outer = observed[
        ["scan_id", "subset", "x_var", "y_var", "n_events", "k_best", "delta_chi2"]
    ].copy()
    observed_outer["outer_global_p"] = observed_outer_p

    if args.output_prefix is None:
        safe_model = re.sub(r"[^A-Za-z0-9_]+", "_", args.null_model)
        prefix = Path(f"tables/gwtc4_population_{safe_model}_{args.outer_n}")
    else:
        prefix = Path(args.output_prefix)

    matrix_path = Path(str(prefix) + "_null_matrix.csv")
    observed_path_out = Path(str(prefix) + "_observed_outercal.csv")
    delta_path = Path(str(prefix) + "_null_delta_matrix.csv")
    metadata_path = Path(str(prefix) + "_metadata.json")
    global_result_path = Path(str(prefix) + "_global_result.json")
    for p in (matrix_path, observed_path_out, delta_path, metadata_path):
        p.parent.mkdir(parents=True, exist_ok=True)

    matrix_df.to_csv(matrix_path, index=False)
    delta_df.to_csv(delta_path, index=False)
    observed_outer.to_csv(observed_path_out, index=False)

    metadata = {
        "schema": "GWTC4_POPULATION_NULL_SENSITIVITY_V1",
        "null_model": args.null_model,
        "null_definition": null_meta,
        "event_summary": str(event_path),
        "observed_manifest": str(args.observed),
        "n_events": int(len(event_df)),
        "n_finite_final_spin": int(len(finite_idx)),
        "n_scans": int(len(observed)),
        "n_unique_subsets": int(len(manifest_subsets)),
        "n_dynamic_final_spin_subsets": int(
            sum(is_dynamic_subset(s) for s in manifest_subsets)
        ),
        "k_min": float(args.k_min),
        "k_max": float(args.k_max),
        "n_k": int(args.n_k),
        "degree": int(args.degree),
        "min_events": int(args.min_events),
        "top_ns": top_ns,
        "seed": int(args.seed),
        "outer_n": int(args.outer_n),
        "observed_reproduction_max_abs_k_error": max_k_err,
        "observed_reproduction_max_abs_delta_error": max_d_err,
        "observed_reproduction_pass": True,
        "observed_outer_p_min": float(np.min(observed_outer_p)),
        "observed_outer_p_median": float(np.median(observed_outer_p)),
        "observed_outer_p_lt_0p05": int(np.sum(observed_outer_p < 0.05)),
        "observed_outer_p_lt_0p10": int(np.sum(observed_outer_p < 0.10)),
        "outer_per_scan_resolution_floor": float(1.0 / (1.0 + args.outer_n)),
        "null_row_rank_resolution_floor": float(1.0 / args.outer_n),
        "output_null_p_matrix": str(matrix_path),
        "output_null_delta_matrix": str(delta_path),
        "output_observed_outercal": str(observed_path_out),
        "suggested_global_result": str(global_result_path),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    print("\nGenerated")
    print("---------")
    print(f"Observed outer-cal vector : {observed_path_out}")
    print(f"Null p-value matrix       : {matrix_path}")
    print(f"Null Delta matrix         : {delta_path}")
    print(f"Metadata                  : {metadata_path}")
    print(f"Per-scan p floor          : {1.0 / (1.0 + args.outer_n):.8g}")
    print(f"Observed p<0.05 count     : {int(np.sum(observed_outer_p < 0.05))}")
    print(f"Observed p<0.10 count     : {int(np.sum(observed_outer_p < 0.10))}")
    print("\nNext command:")
    print(
        "python scripts/run_gwtc_population_global_null.py "
        f"--observed {observed_path_out} --p-column outer_global_p "
        f"--null-matrix {matrix_path} --output {global_result_path}"
    )


if __name__ == "__main__":
    main()
