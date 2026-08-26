#!/usr/bin/env python3
"""Evaluate a frozen GWTC mode on the holdout exactly once, without refitting.

Inputs:
- canonical event table;
- frozen train/holdout manifest;
- frozen model JSON created by ``fit_frozen_gwtc_mode.py``.

This program does NOT scan k, refit the baseline, refit phase/amplitude, or alter
bin edges.  It compares two fully frozen predictive shapes on holdout counts:
(1) smooth baseline and (2) baseline multiplied by the training-selected mode.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_baseline_v2_cv as V2  # noqa: E402


def predictive_delta(counts: np.ndarray, p0: np.ndarray, p1: np.ndarray) -> float:
    counts = np.asarray(counts, float)
    p0 = np.clip(np.asarray(p0, float), 1e-300, None)
    p1 = np.clip(np.asarray(p1, float), 1e-300, None)
    p0 = p0 / p0.sum()
    p1 = p1 / p1.sum()
    # 2*log[L(frozen residual)/L(frozen baseline)] for multinomial counts.
    return float(2.0 * np.sum(counts * np.log(p1 / p0)))


def fixed_model_null_p(
    counts: np.ndarray,
    p0: np.ndarray,
    p1: np.ndarray,
    null_n: int,
    seed: int,
) -> tuple[float, np.ndarray]:
    n = int(np.sum(counts))
    if n <= 0:
        raise ValueError("Holdout contains zero in-support events")
    rng = np.random.default_rng(seed)
    p0 = np.clip(np.asarray(p0, float), 1e-300, None)
    p0 = p0 / p0.sum()
    observed = predictive_delta(counts, p0, p1)
    null = np.empty(int(null_n), float)
    for i in range(int(null_n)):
        y = rng.multinomial(n, p0)
        null[i] = predictive_delta(y, p0, p1)
    p = (1 + int(np.sum(null >= observed))) / (1 + int(null_n))
    return float(p), null


def extract_holdout_values(
    df: pd.DataFrame,
    manifest: pd.DataFrame,
    variable: str,
) -> np.ndarray:
    holdout = V2.apply_manifest_split(df, manifest, "holdout")
    if variable not in holdout.columns:
        raise ValueError(f"Variable {variable!r} not present in event table")
    values = pd.to_numeric(holdout[variable], errors="coerce").to_numpy(float)
    return values[np.isfinite(values) & (values > 0)]


def evaluate(
    values: np.ndarray,
    model: dict,
    null_n: int,
    seed: int,
) -> dict:
    values = np.asarray(values, float)
    edges = np.asarray(model["bin_edges"], float)
    p0 = np.asarray(model["baseline_probabilities"], float)
    p1 = np.asarray(model["residual_probabilities"], float)
    if len(edges) != len(p0) + 1 or len(p0) != len(p1):
        raise ValueError("Frozen model has inconsistent bin/probability dimensions")

    ell = np.log(values)
    lo = float(edges[0])
    hi = float(edges[-1])
    in_support = (ell >= lo) & (ell <= hi)
    outside = int(np.sum(~in_support))
    used = ell[in_support]
    counts, _ = np.histogram(used, bins=edges)

    delta = predictive_delta(counts, p0, p1)
    p_value, null = fixed_model_null_p(counts, p0, p1, null_n=null_n, seed=seed)

    if outside > 0:
        verdict = "INCOMPLETE_OUT_OF_SUPPORT"
    elif delta <= 0:
        verdict = "FAIL"
    elif p_value < 0.05:
        verdict = "PASS_FIXED_MODE"
    else:
        verdict = "PARTIAL"

    return {
        "protocol": "GWTC_STRICT_FROZEN_HOLDOUT_V1",
        "variable": model["variable"],
        "baseline_degree": int(model["baseline_degree"]),
        "bins": int(model["bins"]),
        "training_n": int(model["training_n"]),
        "holdout_positive_n": int(len(values)),
        "holdout_in_support_n": int(len(used)),
        "holdout_outside_support_n": outside,
        "k_frozen": float(model["k_star"]),
        "n_frozen_training_support": float(model["n_star"]),
        "phase_frozen": float(model["residual_phase"]),
        "amplitude_frozen": float(model["residual_amplitude"]),
        "holdout_delta_deviance": float(delta),
        "holdout_log_likelihood_ratio": float(delta / 2.0),
        "null_n": int(null_n),
        "null_seed": int(seed),
        "holdout_p": float(p_value),
        "null_mean_delta": float(np.mean(null)),
        "null_sd_delta": float(np.std(null, ddof=1)) if len(null) > 1 else 0.0,
        "verdict": verdict,
        "class_i_claim": False,
        "notes": (
            "PASS_FIXED_MODE, if reached, means only that the frozen training-selected "
            "mode predicts the holdout better than the frozen smooth baseline under this "
            "test. It is not a WCT Class-I claim."
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--table", default="tables/gwtc_events_clean.csv")
    p.add_argument("--manifest", default="tables/gwtc_holdout_manifest.csv")
    p.add_argument("--model", default="outputs/summary/gwtc_frozen_mode.json")
    p.add_argument("--null-n", type=int, default=10000)
    p.add_argument("--seed", type=int, default=271828)
    p.add_argument("--output", default="tables/gwtc_frozen_holdout_result.csv")
    args = p.parse_args()

    model = json.loads(Path(args.model).read_text(encoding="utf-8"))
    if model.get("protocol") != "GWTC_STRICT_FROZEN_MODE_V1":
        raise SystemExit("Refusing non-strict or unknown frozen model protocol")
    if model.get("holdout_evaluated") is not False:
        raise SystemExit("Frozen model is not marked holdout_evaluated=false")

    df = pd.read_csv(args.table)
    manifest = pd.read_csv(args.manifest)
    try:
        values = extract_holdout_values(df, manifest, model["variable"])
        result = evaluate(values, model, null_n=args.null_n, seed=args.seed)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([result]).to_csv(out, index=False)

    print(f"Wrote {out}")
    print(f"Frozen k       : {result['k_frozen']:.8g}")
    print(f"Frozen phase   : {result['phase_frozen']:.8g}")
    print(f"Frozen amplitude: {result['amplitude_frozen']:.8g}")
    print(f"Holdout N      : {result['holdout_in_support_n']} in support; "
          f"{result['holdout_outside_support_n']} outside")
    print(f"DeltaD holdout : {result['holdout_delta_deviance']:.8g}")
    print(f"p holdout      : {result['holdout_p']:.8g} (N_null={result['null_n']})")
    print(f"VERDICT        : {result['verdict']}")


if __name__ == "__main__":
    main()
