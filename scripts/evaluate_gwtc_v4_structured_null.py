#!/usr/bin/env python3
"""Evaluate the frozen GWTC V4 fixed-k structured-null robustness challenge.

This evaluator reads a previously frozen V4 JSON.  It cannot fit a population,
change k, refit residual coefficients, or select a new baseline.  It computes
one observed GWTC-5 fixed-model statistic at the externally published k=9.7
and compares it against Monte Carlo catalogs drawn from each frozen structured
non-periodic null scenario.

Because GWTC-5 was already inspected in V2/V3, V4 is a robustness challenge,
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
from gwtc_structured_null import (  # noqa: E402
    StructuredNullParams,
    structured_fixed_statistic_null,
)
from gwtc_unbinned_kde import fixed_log_likelihood_ratio  # noqa: E402

EXPECTED_K = 9.7


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
        raise SystemExit("Clean table lacks is_primary_entry; refusing noncanonical V4 evaluation.")
    s = df["is_primary_entry"]
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--table", default="tables/gwtc_events_clean.csv")
    p.add_argument("--manifest", default="tables/gwtc_holdout_manifest.csv")
    p.add_argument("--frozen", default="tables/gwtc_v4_frozen_structured_null.json")
    p.add_argument("--null-n", type=int, default=10000)
    p.add_argument("--seed", type=int, default=2718281)
    p.add_argument("--grid-n", type=int, default=4096)
    p.add_argument("--output", default="tables/gwtc_v4_structured_null_result.csv")
    args = p.parse_args()

    if args.null_n <= 0:
        raise SystemExit("--null-n must be positive")

    table_path = Path(args.table)
    manifest_path = Path(args.manifest)
    frozen_path = Path(args.frozen)
    frozen = json.loads(frozen_path.read_text(encoding="utf-8"))

    if frozen.get("schema") != "GWTC_V4_STRUCTURED_NULL_FIXED_PREDICTION_V1":
        raise SystemExit("Unrecognized frozen V4 schema.")
    if frozen.get("prediction_status") != "EXTERNALLY_PUBLISHED_AND_FROZEN_BEFORE_V4":
        raise SystemExit("V4 prediction-status marker missing; refusing evaluation.")
    if abs(float(frozen.get("k_pred", np.nan)) - EXPECTED_K) > 1e-12:
        raise SystemExit("Frozen V4 k is not the published k=9.7 prediction.")
    if sha256_file(manifest_path) != frozen.get("manifest_sha256"):
        raise SystemExit("Manifest hash differs from frozen V4; refusing evaluation.")
    if sha256_file(table_path) != frozen.get("clean_table_sha256"):
        raise SystemExit("Clean-table hash differs from frozen V4; refusing evaluation.")

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

    mode = frozen["fixed_mode"]
    k = float(mode["k"])
    a = float(mode["a_cos"])
    b = float(mode["b_sin"])
    normalizer = float(mode["normalizer_Z"])
    if abs(k - EXPECTED_K) > 1e-12:
        raise SystemExit("Fixed-mode k differs from the published k=9.7.")

    observed = fixed_log_likelihood_ratio(ell, k, a, b, normalizer)

    scenarios = frozen.get("structured_null_scenarios", [])
    if not scenarios:
        raise SystemExit("Frozen V4 contains no structured-null scenarios.")

    rows = []
    p_values = []
    for index, scenario in enumerate(scenarios):
        params = StructuredNullParams.from_dict(scenario["params"])
        scenario_seed = int(args.seed) + index * 1000003
        null = structured_fixed_statistic_null(
            params=params,
            n_events=int(ell.size),
            k=k,
            a=a,
            b=b,
            normalizer=normalizer,
            null_n=int(args.null_n),
            seed=scenario_seed,
            grid_n=int(args.grid_n),
        )
        n_ge = int(np.sum(null >= observed))
        p_value = float((1 + n_ge) / (1 + int(args.null_n)))
        p_values.append(p_value)

        if args.null_n < 1000:
            verdict = "PARTIAL_NULL_TOO_SMALL"
        elif observed > 0.0 and p_value < 0.05:
            verdict = "PASS_STRUCTURED_NULL_SCENARIO"
        else:
            verdict = "FAIL_STRUCTURED_NULL_SCENARIO"

        rows.append(
            {
                "schema": "GWTC_V4_STRUCTURED_NULL_RESULT_V1",
                "scientific_status": "ROBUSTNESS_CHALLENGE_NOT_INDEPENDENT_REPLICATION",
                "variable": variable,
                "test": "fixed_k_9p7_unbinned_residual_statistic_vs_structured_nonperiodic_null",
                "frozen_v4_sha256": sha256_file(frozen_path),
                "manifest_sha256": sha256_file(manifest_path),
                "clean_table_sha256": sha256_file(table_path),
                "holdout_manifest_rows": int((manifest["split"] == "holdout").sum()),
                "holdout_n_positive_finite": int(ell.size),
                "scenario": str(scenario["name"]),
                "population_family": str(scenario["population_family"]),
                "selection_family": str(scenario["selection_family"]),
                "selection_gamma": float(scenario["selection_gamma"]),
                "structured_null_params_json": json.dumps(scenario["params"], sort_keys=True, separators=(",", ":")),
                "k_fixed": k,
                "a_cos_frozen": a,
                "b_sin_frozen": b,
                "amplitude_frozen": float(mode["amplitude"]),
                "normalizer_Z_frozen": normalizer,
                "delta_2logl_holdout": observed,
                "null_n": int(args.null_n),
                "seed": scenario_seed,
                "null_ge_observed": n_ge,
                "p_structured_null": p_value,
                "verdict": verdict,
                "interpretation_limit": "Phenomenological broken-power-law+peak selection challenge; not full hierarchical LVK population inference or an independent future-catalog replication.",
            }
        )

    worst_p = float(max(p_values))
    if args.null_n < 1000:
        overall = "PARTIAL_NULL_TOO_SMALL"
    elif observed > 0.0 and worst_p < 0.05:
        overall = "PASS_ALL_STRUCTURED_NULLS"
    else:
        overall = "FAIL_AT_LEAST_ONE_STRUCTURED_NULL"

    for row in rows:
        row["worst_case_p_across_scenarios"] = worst_p
        row["overall_verdict"] = overall

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)

    print(f"Wrote {out}")
    print(f"Published k      : {k:.8g} (NO SCAN)")
    print(f"Frozen amplitude : {float(mode['amplitude']):.8g}")
    print(f"Holdout N        : {ell.size} positive finite")
    print(f"Delta 2logL      : {observed:.8g}")
    for row in rows:
        print(
            f"gamma={row['selection_gamma']:g}: "
            f"p={row['p_structured_null']:.8g}, "
            f"ge={row['null_ge_observed']}/{row['null_n']}, "
            f"{row['verdict']}"
        )
    print(f"Worst-case p     : {worst_p:.8g}")
    print(f"VERDICT          : {overall}")


if __name__ == "__main__":
    main()
