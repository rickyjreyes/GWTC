#!/usr/bin/env python3
"""Fit and freeze the GWTC V3 unbinned KDE residual model on TRAIN only.

This script must be run before the V3 holdout evaluator.  It reads the already
frozen V2-style manifest, uses only rows assigned to split=train, selects a
non-polynomial Gaussian-KDE bandwidth by leave-one-out likelihood, scans the
declared k grid on training data, and writes every parameter needed for a later
no-refit holdout test to a tracked JSON artifact.

GWTC-5 has already been inspected by the earlier V2 experiment, so V3 is a
predeclared robustness test, not a new independent catalog replication.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gwtc_unbinned_kde import scan_training_mode, select_bandwidth_loo  # noqa: E402


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def event_column(df: pd.DataFrame) -> str:
    for name in ("commonName", "common_name", "event", "event_name", "name"):
        if name in df.columns:
            return name
    raise SystemExit("No event-name column found in clean table.")


def primary_mask(df: pd.DataFrame) -> pd.Series:
    if "is_primary_entry" not in df.columns:
        raise SystemExit("Clean table lacks is_primary_entry; refusing noncanonical V3 fit.")
    s = df["is_primary_entry"]
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--table", default="tables/gwtc_events_clean.csv")
    p.add_argument("--manifest", default="tables/gwtc_holdout_manifest.csv")
    p.add_argument("--variable", default="M_chirp")
    p.add_argument("--bandwidth-multipliers", default="0.5,0.75,1.0,1.25,1.5,2.0,2.5")
    p.add_argument("--k-min", type=float, default=0.5)
    p.add_argument("--k-max", type=float, default=40.0)
    p.add_argument("--n-k", type=int, default=120)
    p.add_argument("--gh-n", type=int, default=40)
    p.add_argument("--output", default="tables/gwtc_v3_frozen_unbinned_kde_mode.json")
    args = p.parse_args()

    table_path = Path(args.table)
    manifest_path = Path(args.manifest)
    df = pd.read_csv(table_path)
    manifest = pd.read_csv(manifest_path)

    if args.variable not in df.columns:
        raise SystemExit(f"Variable {args.variable!r} not found in {args.table}")
    if not {"event_name", "split"}.issubset(manifest.columns):
        raise SystemExit("Manifest must contain event_name and split columns.")
    if not (manifest["split"] == "train").any() or not (manifest["split"] == "holdout").any():
        raise SystemExit("Manifest must already contain both train and holdout rows.")

    train_names = set(manifest.loc[manifest["split"] == "train", "event_name"].astype(str))
    ev_col = event_column(df)
    canonical = df.loc[primary_mask(df)].copy()
    canonical[ev_col] = canonical[ev_col].astype(str)
    train = canonical[canonical[ev_col].isin(train_names)].copy()

    z = pd.to_numeric(train[args.variable], errors="coerce").to_numpy(float)
    z = z[np.isfinite(z) & (z > 0.0)]
    if z.size < 30:
        raise SystemExit(f"Too few positive finite training {args.variable} values: {z.size}")
    ell = np.log(z)

    multipliers = [float(x.strip()) for x in args.bandwidth_multipliers.split(",") if x.strip()]
    bw = select_bandwidth_loo(ell, multipliers)

    if args.n_k < 2 or args.k_max <= args.k_min:
        raise SystemExit("Invalid k grid.")
    k_grid = np.linspace(args.k_min, args.k_max, args.n_k)
    best, scan = scan_training_mode(
        ell,
        ell,  # KDE centers are exactly the frozen training log-values.
        bw.selected_bandwidth,
        k_grid,
        gh_n=args.gh_n,
    )

    delta_ell = float(np.max(ell) - np.min(ell))
    n_equiv = float(best.k * delta_ell / (2.0 * math.pi))

    payload = {
        "schema": "GWTC_V3_UNBINNED_KDE_FIXED_MODE_V1",
        "scientific_status": "ROBUSTNESS_TEST_NOT_INDEPENDENT_REPLICATION",
        "reason": "GWTC-5 holdout was previously inspected by V2; V3 tests survival under a new predeclared non-polynomial unbinned analysis.",
        "variable": args.variable,
        "coordinate": "ell=log(variable)",
        "baseline_family": "gaussian_kde",
        "bandwidth_selection": "training-only leave-one-out log likelihood",
        "manifest_sha256": sha256_file(manifest_path),
        "clean_table_sha256": sha256_file(table_path),
        "manifest_path": str(manifest_path).replace("\\", "/"),
        "clean_table_path": str(table_path).replace("\\", "/"),
        "training_manifest_rows": int((manifest["split"] == "train").sum()),
        "training_n_positive_finite": int(ell.size),
        "train_ell": [float(x) for x in ell],
        "train_ell_min": float(np.min(ell)),
        "train_ell_max": float(np.max(ell)),
        "train_delta_ell": delta_ell,
        "scott_bandwidth": float(bw.base_scott),
        "bandwidth_candidates": [float(x) for x in multipliers],
        "bandwidth_loo_loglik": {str(k): float(v) for k, v in bw.scores.items()},
        "selected_bandwidth_multiplier": float(bw.selected_multiplier),
        "selected_bandwidth": float(bw.selected_bandwidth),
        "k_min": float(args.k_min),
        "k_max": float(args.k_max),
        "n_k": int(args.n_k),
        "gh_n": int(args.gh_n),
        "coefficient_bounds": [-3.0, 3.0],
        "k_star": float(best.k),
        "a_cos": float(best.a),
        "b_sin": float(best.b),
        "amplitude": float(best.amplitude),
        "phase_atan2_b_a": float(best.phase_atan2_b_a),
        "normalizer_Z": float(best.normalizer),
        "n_equiv_over_training_span": n_equiv,
        "train_delta_2logl": float(best.train_delta_2logl),
        "scan": [
            {
                "k": float(row.k),
                "a_cos": float(row.a),
                "b_sin": float(row.b),
                "amplitude": float(row.amplitude),
                "phase_atan2_b_a": float(row.phase_atan2_b_a),
                "normalizer_Z": float(row.normalizer),
                "train_delta_2logl": float(row.train_delta_2logl),
            }
            for row in scan
        ],
        "holdout_evaluated": False,
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"Wrote {out}")
    print(f"Training N       : {ell.size}")
    print(f"KDE bandwidth    : {bw.selected_bandwidth:.8g} (multiplier={bw.selected_multiplier:g})")
    print(f"Frozen k         : {best.k:.8g}")
    print(f"Frozen a,b       : {best.a:.8g}, {best.b:.8g}")
    print(f"Frozen amplitude : {best.amplitude:.8g}")
    print(f"Equivalent n     : {n_equiv:.8g}")
    print(f"Train Delta2logL : {best.train_delta_2logl:.8g}")
    print("HOLDOUT HAS NOT BEEN EVALUATED BY THIS SCRIPT.")
    print("Hash/commit this JSON before running the V3 evaluator.")


if __name__ == "__main__":
    main()
