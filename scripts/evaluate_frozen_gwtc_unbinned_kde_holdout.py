#!/usr/bin/env python3
"""Evaluate the frozen GWTC V3 unbinned KDE mode on the declared holdout.

This program cannot scan k, select a bandwidth, or refit a/b.  It reads those
objects from the frozen training JSON and computes a single fixed-model
unbinned likelihood-ratio statistic on the manifest holdout.  Its null is
sampled from the exact frozen Gaussian-KDE mixture and nothing is refit in a
null replicate.

Because GWTC-5 was already inspected by V2, a PASS here is a robustness PASS,
not an independent future-catalog replication.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gwtc_unbinned_kde import fixed_log_likelihood_ratio, fixed_model_null  # noqa: E402


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
        raise SystemExit("Clean table lacks is_primary_entry; refusing noncanonical V3 evaluation.")
    s = df["is_primary_entry"]
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--table", default="tables/gwtc_events_clean.csv")
    p.add_argument("--manifest", default="tables/gwtc_holdout_manifest.csv")
    p.add_argument("--frozen", default="tables/gwtc_v3_frozen_unbinned_kde_mode.json")
    p.add_argument("--null-n", type=int, default=10000)
    p.add_argument("--seed", type=int, default=314159)
    p.add_argument("--output", default="tables/gwtc_v3_unbinned_kde_holdout_result.csv")
    args = p.parse_args()

    table_path = Path(args.table)
    manifest_path = Path(args.manifest)
    frozen_path = Path(args.frozen)
    frozen = json.loads(frozen_path.read_text(encoding="utf-8"))

    if frozen.get("schema") != "GWTC_V3_UNBINNED_KDE_FIXED_MODE_V1":
        raise SystemExit("Unrecognized frozen V3 schema.")
    if sha256_file(manifest_path) != frozen.get("manifest_sha256"):
        raise SystemExit("Manifest hash differs from the one used to freeze V3; refusing evaluation.")
    if sha256_file(table_path) != frozen.get("clean_table_sha256"):
        raise SystemExit("Clean-table hash differs from the one used to freeze V3; refusing evaluation.")

    df = pd.read_csv(table_path)
    manifest = pd.read_csv(manifest_path)
    if not {"event_name", "split"}.issubset(manifest.columns):
        raise SystemExit("Manifest must contain event_name and split columns.")

    variable = str(frozen["variable"])
    if variable not in df.columns:
        raise SystemExit(f"Variable {variable!r} not found in clean table.")

    holdout_names = set(manifest.loc[manifest["split"] == "holdout", "event_name"].astype(str))
    ev_col = event_column(df)
    canonical = df.loc[primary_mask(df)].copy()
    canonical[ev_col] = canonical[ev_col].astype(str)
    holdout = canonical[canonical[ev_col].isin(holdout_names)].copy()

    z = pd.to_numeric(holdout[variable], errors="coerce").to_numpy(float)
    z = z[np.isfinite(z) & (z > 0.0)]
    if z.size == 0:
        raise SystemExit("No positive finite holdout values.")
    ell = np.log(z)

    centers = np.asarray(frozen["train_ell"], dtype=float)
    bandwidth = float(frozen["selected_bandwidth"])
    k = float(frozen["k_star"])
    a = float(frozen["a_cos"])
    b = float(frozen["b_sin"])
    normalizer = float(frozen["normalizer_Z"])

    observed = fixed_log_likelihood_ratio(ell, k, a, b, normalizer)
    null = fixed_model_null(
        kde_centers=centers,
        bandwidth=bandwidth,
        n_events=int(ell.size),
        k=k,
        a=a,
        b=b,
        normalizer=normalizer,
        null_n=args.null_n,
        seed=args.seed,
    )
    n_ge = int(np.sum(null >= observed))
    p_value = float((1 + n_ge) / (1 + args.null_n))

    lo = float(frozen["train_ell_min"])
    hi = float(frozen["train_ell_max"])
    n_inside = int(np.sum((ell >= lo) & (ell <= hi)))
    n_outside = int(ell.size - n_inside)

    if args.null_n < 1000:
        verdict = "PARTIAL_NULL_TOO_SMALL"
    elif observed > 0.0 and p_value < 0.05:
        verdict = "PASS_ROBUSTNESS_FIXED_MODE"
    else:
        verdict = "FAIL_ROBUSTNESS_FIXED_MODE"

    row = {
        "schema": "GWTC_V3_UNBINNED_KDE_HOLDOUT_RESULT_V1",
        "scientific_status": "ROBUSTNESS_TEST_NOT_INDEPENDENT_REPLICATION",
        "variable": variable,
        "baseline_family": "gaussian_kde_loo_cv",
        "test": "unbinned_fixed_model_2log_likelihood_ratio",
        "frozen_model_sha256": sha256_file(frozen_path),
        "manifest_sha256": sha256_file(manifest_path),
        "clean_table_sha256": sha256_file(table_path),
        "holdout_manifest_rows": int((manifest["split"] == "holdout").sum()),
        "holdout_n_positive_finite": int(ell.size),
        "holdout_n_inside_training_span": n_inside,
        "holdout_n_outside_training_span": n_outside,
        "selected_bandwidth": bandwidth,
        "k_frozen": k,
        "a_cos_frozen": a,
        "b_sin_frozen": b,
        "amplitude_frozen": float(frozen["amplitude"]),
        "normalizer_Z_frozen": normalizer,
        "delta_2logl_holdout": observed,
        "null_n": int(args.null_n),
        "seed": int(args.seed),
        "null_ge_observed": n_ge,
        "p_holdout": p_value,
        "verdict": verdict,
        "interpretation_limit": "V3 is a robustness test because GWTC-5 was already inspected in V2; it is not an independent future-catalog replication.",
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(out, index=False)

    print(f"Wrote {out}")
    print(f"Frozen KDE h     : {bandwidth:.8g}")
    print(f"Frozen k         : {k:.8g}")
    print(f"Frozen amplitude : {float(frozen['amplitude']):.8g}")
    print(f"Holdout N        : {ell.size} positive finite")
    print(f"Training-span    : {n_inside} inside; {n_outside} outside")
    print(f"Delta 2logL      : {observed:.8g}")
    print(f"p holdout        : {p_value:.8g} (N_null={args.null_n})")
    print(f"VERDICT          : {verdict}")


if __name__ == "__main__":
    main()
