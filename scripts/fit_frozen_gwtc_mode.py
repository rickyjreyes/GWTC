#!/usr/bin/env python3
"""Fit and freeze one GWTC residual mode using the TRAINING split only.

This program never evaluates holdout residuals.  It reads a pre-existing frozen
manifest, selects only ``split=train`` events, reads the already cross-validated
baseline degree, fits the smooth baseline, scans the declared resolvable k grid,
and writes the resulting baseline + residual prediction to JSON.

The JSON is the immutable object consumed by ``evaluate_frozen_gwtc_holdout``.
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
import wct_gwtc_lib as L  # noqa: E402


def selected_degree(path: str) -> int:
    cv = pd.read_csv(path)
    if "selected" not in cv.columns or "degree" not in cv.columns:
        raise ValueError("Baseline CV file must contain selected and degree columns")
    selected = cv[cv["selected"].astype(str).str.lower().isin(("true", "1"))]
    if "data_split" in selected.columns:
        selected = selected[selected["data_split"].astype(str) == "train"]
    if len(selected) != 1:
        raise ValueError(f"Expected exactly one selected TRAIN baseline row, found {len(selected)}")
    return int(selected.iloc[0]["degree"])


def fit_training_mode(
    values: np.ndarray,
    baseline_degree: int,
    bins: int,
    k_min: float,
    k_max: float,
    n_k: int,
) -> dict:
    values = np.asarray(values, float)
    values = values[np.isfinite(values) & (values > 0)]
    if len(values) < max(30, bins):
        raise ValueError(f"Too few positive training values: {len(values)}")

    ell = np.log(values)
    lo = float(np.min(ell))
    hi = float(np.max(ell))
    if hi <= lo:
        raise ValueError("Degenerate training support")

    edges = np.linspace(lo, hi, bins + 1)
    counts, _ = np.histogram(ell, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    center = float(np.mean(centers))
    scale = float(np.std(centers))
    if scale <= 0 or not np.isfinite(scale):
        raise ValueError("Degenerate baseline scale")

    X0 = V2.design_poly(centers, baseline_degree, center, scale)
    beta0 = V2.fit_poisson_glm(X0, counts)
    mu0 = np.exp(np.clip(X0 @ beta0, -30.0, 30.0))

    delta_ell = hi - lo
    nyquist_k = math.pi * bins / delta_ell
    effective_k_max = min(float(k_max), float(nyquist_k))
    if effective_k_max <= k_min:
        raise ValueError("No resolvable k range below the histogram Nyquist frequency")
    k_grid = np.linspace(float(k_min), effective_k_max, int(n_k))

    d0 = V2.poisson_deviance(counts, mu0)
    dd = []
    coeffs = []
    models = []
    offset = np.log(np.clip(mu0, 1e-12, None))
    for k in k_grid:
        Xr = np.column_stack(
            [np.ones_like(centers), np.cos(k * centers), np.sin(k * centers)]
        )
        beta_r, mu_r = L.fit_poisson_glm(counts, Xr, offset=offset)
        dd.append(d0 - V2.poisson_deviance(counts, mu_r))
        coeffs.append(beta_r)
        models.append(mu_r)

    idx = int(np.argmax(dd))
    k_star = float(k_grid[idx])
    delta_d = float(dd[idx])
    c, a, b = [float(x) for x in coeffs[idx]]
    mu1 = np.asarray(models[idx], float)

    # A*cos(k*ell + phi) = a*cos(k*ell) + b*sin(k*ell)
    amplitude = float(math.hypot(a, b))
    phase = float(math.atan2(-b, a))
    n_star = float(k_star * delta_ell / (2.0 * math.pi))

    p0 = np.clip(mu0, 1e-300, None)
    p0 = p0 / p0.sum()
    p1 = np.clip(mu1, 1e-300, None)
    p1 = p1 / p1.sum()

    return {
        "protocol": "GWTC_STRICT_FROZEN_MODE_V1",
        "variable": None,
        "baseline_family": "polynomial_poisson",
        "baseline_degree": int(baseline_degree),
        "bins": int(bins),
        "training_n": int(len(values)),
        "ell_min": lo,
        "ell_max": hi,
        "delta_ell_active": delta_ell,
        "bin_edges": edges.tolist(),
        "bin_centers": centers.tolist(),
        "baseline_center": center,
        "baseline_scale": scale,
        "baseline_beta": beta0.tolist(),
        "baseline_probabilities": p0.tolist(),
        "k_grid_declared_min": float(k_min),
        "k_grid_declared_max": float(k_max),
        "k_grid_effective_max": effective_k_max,
        "nyquist_k": float(nyquist_k),
        "n_k": int(n_k),
        "k_star": k_star,
        "n_star": n_star,
        "residual_c": c,
        "residual_a": a,
        "residual_b": b,
        "residual_amplitude": amplitude,
        "residual_phase": phase,
        "residual_probabilities": p1.tolist(),
        "training_delta_deviance": delta_d,
        "holdout_evaluated": False,
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--table", default="tables/gwtc_events_clean.csv")
    p.add_argument("--manifest", default="tables/gwtc_holdout_manifest.csv")
    p.add_argument("--baseline-cv", default="tables/gwtc_baseline_v2_cv.csv")
    p.add_argument("--variable", default="M_chirp")
    p.add_argument("--bins", type=int, default=40)
    p.add_argument("--k-min", type=float, default=0.5)
    p.add_argument("--k-max", type=float, default=40.0)
    p.add_argument("--n-k", type=int, default=120)
    p.add_argument("--output", default="outputs/summary/gwtc_frozen_mode.json")
    args = p.parse_args()

    df = pd.read_csv(args.table)
    manifest = pd.read_csv(args.manifest)
    if args.variable not in df.columns:
        raise SystemExit(f"Variable {args.variable!r} not found")

    try:
        train = V2.apply_manifest_split(df, manifest, "train")
        degree = selected_degree(args.baseline_cv)
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    values = pd.to_numeric(train[args.variable], errors="coerce").to_numpy(float)
    try:
        model = fit_training_mode(
            values,
            baseline_degree=degree,
            bins=args.bins,
            k_min=args.k_min,
            k_max=args.k_max,
            n_k=args.n_k,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    model["variable"] = args.variable
    model["manifest"] = args.manifest
    model["baseline_cv"] = args.baseline_cv

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(model, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(
        f"Wrote {out}; training_n={model['training_n']}; degree={degree}; "
        f"k*={model['k_star']:.6g}; n*={model['n_star']:.6g}; "
        f"phase={model['residual_phase']:.6g}; A={model['residual_amplitude']:.6g}; "
        f"DeltaD_train={model['training_delta_deviance']:.6g}"
    )
    print("HOLDOUT HAS NOT BEEN EVALUATED. Hash/archive this JSON before the next step.")


if __name__ == "__main__":
    main()
