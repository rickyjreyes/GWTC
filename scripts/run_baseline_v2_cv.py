#!/usr/bin/env python3
"""Cross-validate flexible smooth baselines for GWTC log-domain counts.

This is deliberately a baseline-only challenge. It does not scan or fit a WCT
residual. Hyperparameters are selected from training folds using Poisson
predictive deviance so residual testing can occur only after the smooth
background model is frozen.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize


def poisson_deviance(y: np.ndarray, mu: np.ndarray) -> float:
    mu = np.clip(np.asarray(mu, float), 1e-12, None)
    y = np.asarray(y, float)
    term = np.zeros_like(y)
    positive = y > 0
    term[positive] = y[positive] * np.log(y[positive] / mu[positive])
    return float(2.0 * np.sum(term - (y - mu)))


def design_poly(x: np.ndarray, degree: int, center: float, scale: float) -> np.ndarray:
    t = (x - center) / scale
    return np.column_stack([t ** d for d in range(degree + 1)])


def fit_poisson_glm(X: np.ndarray, y: np.ndarray) -> np.ndarray:
    y = np.asarray(y, float)

    def objective(beta: np.ndarray) -> float:
        eta = np.clip(X @ beta, -30.0, 30.0)
        mu = np.exp(eta)
        return float(np.sum(mu - y * eta))

    beta0 = np.zeros(X.shape[1], dtype=float)
    beta0[0] = math.log(max(float(np.mean(y)), 1e-6))
    result = minimize(objective, beta0, method="BFGS")
    if not result.success and not np.isfinite(result.fun):
        raise RuntimeError(f"Poisson GLM fit failed: {result.message}")
    return np.asarray(result.x, float)


def make_folds(n: int, n_folds: int) -> list[np.ndarray]:
    # Interleaved deterministic folds preserve coverage across the log domain.
    indices = np.arange(n)
    return [indices[fold::n_folds] for fold in range(n_folds)]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", default="tables/gwtc_events_clean.csv")
    p.add_argument("--variable", default="M_chirp")
    p.add_argument("--bins", type=int, default=40)
    p.add_argument("--degrees", default="3,4,5,6,7,8")
    p.add_argument("--folds", type=int, default=5)
    p.add_argument("--p-astro-min", type=float, default=0.5)
    p.add_argument("--output", default="tables/gwtc_baseline_v2_cv.csv")
    args = p.parse_args()

    df = pd.read_csv(args.input)
    if args.variable not in df.columns:
        raise SystemExit(f"Variable {args.variable!r} not found in {args.input}")
    if "p_astro" in df.columns:
        df = df[pd.to_numeric(df["p_astro"], errors="coerce") >= args.p_astro_min]

    z = pd.to_numeric(df[args.variable], errors="coerce").to_numpy(float)
    z = z[np.isfinite(z) & (z > 0)]
    if len(z) < max(30, args.bins):
        raise SystemExit(f"Too few positive finite {args.variable} values: {len(z)}")

    ell = np.log(z)
    counts, edges = np.histogram(ell, bins=args.bins)
    centers = 0.5 * (edges[:-1] + edges[1:])
    center = float(np.mean(centers))
    scale = float(np.std(centers))
    if not np.isfinite(scale) or scale <= 0:
        raise SystemExit("Degenerate log-domain support")

    degrees = [int(x.strip()) for x in args.degrees.split(",") if x.strip()]
    folds = make_folds(len(counts), args.folds)
    rows: list[dict[str, float | int | str]] = []

    for degree in degrees:
        fold_scores = []
        for fold_idx, valid_idx in enumerate(folds):
            train_mask = np.ones(len(counts), dtype=bool)
            train_mask[valid_idx] = False
            train_idx = np.flatnonzero(train_mask)

            X_train = design_poly(centers[train_idx], degree, center, scale)
            beta = fit_poisson_glm(X_train, counts[train_idx])
            X_valid = design_poly(centers[valid_idx], degree, center, scale)
            mu_valid = np.exp(np.clip(X_valid @ beta, -30.0, 30.0))
            dev = poisson_deviance(counts[valid_idx], mu_valid)
            fold_scores.append(dev)
            rows.append(
                {
                    "variable": args.variable,
                    "bins": args.bins,
                    "baseline_family": "polynomial_poisson",
                    "degree": degree,
                    "fold": fold_idx,
                    "validation_deviance": dev,
                    "n_events": len(z),
                    "p_astro_min": args.p_astro_min,
                }
            )

        rows.append(
            {
                "variable": args.variable,
                "bins": args.bins,
                "baseline_family": "polynomial_poisson",
                "degree": degree,
                "fold": "mean",
                "validation_deviance": float(np.mean(fold_scores)),
                "n_events": len(z),
                "p_astro_min": args.p_astro_min,
            }
        )

    out = pd.DataFrame(rows)
    means = out[out["fold"].astype(str) == "mean"].copy()
    winner = means.loc[means["validation_deviance"].astype(float).idxmin()]
    out["selected"] = False
    out.loc[
        (out["degree"] == int(winner["degree"])) & (out["fold"].astype(str) == "mean"),
        "selected",
    ] = True

    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    print(
        f"Wrote {path}; selected polynomial degree={int(winner['degree'])} "
        f"by mean validation deviance={float(winner['validation_deviance']):.6g}"
    )


if __name__ == "__main__":
    main()
