#!/usr/bin/env python3
"""Deterministic helpers for the GWTC V3 unbinned KDE robustness test.

The baseline is a one-dimensional Gaussian kernel density estimate in
ell = log(z).  Its bandwidth is selected on training data only by leave-one-out
log likelihood.  A residual mode is an exponential tilt of that fixed density:

    f1(ell) = f0(ell) * exp(a cos(k ell) + b sin(k ell)) / Z.

This module contains no catalog-selection logic and never inspects a holdout.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize

SQRT_2PI = math.sqrt(2.0 * math.pi)


def scott_bandwidth(x: Sequence[float]) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        raise ValueError("Need at least three finite points for a KDE bandwidth.")
    scale = float(np.std(x, ddof=1))
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("Degenerate sample scale for KDE bandwidth.")
    return scale * x.size ** (-1.0 / 5.0)


def loo_loglik_gaussian_kde(x: Sequence[float], bandwidth: float) -> float:
    """Exact O(n^2) leave-one-out Gaussian-KDE log likelihood."""
    x = np.asarray(x, dtype=float)
    if bandwidth <= 0.0 or not np.isfinite(bandwidth):
        raise ValueError("bandwidth must be positive and finite")
    n = x.size
    if n < 3:
        raise ValueError("Need at least three points for leave-one-out KDE.")
    u = (x[:, None] - x[None, :]) / bandwidth
    kernel = np.exp(-0.5 * u * u) / (SQRT_2PI * bandwidth)
    np.fill_diagonal(kernel, 0.0)
    density = kernel.sum(axis=1) / float(n - 1)
    density = np.clip(density, 1e-300, None)
    return float(np.sum(np.log(density)))


@dataclass(frozen=True)
class BandwidthSelection:
    base_scott: float
    selected_multiplier: float
    selected_bandwidth: float
    scores: dict[float, float]


def select_bandwidth_loo(
    x: Sequence[float], multipliers: Sequence[float]
) -> BandwidthSelection:
    x = np.asarray(x, dtype=float)
    base = scott_bandwidth(x)
    scores: dict[float, float] = {}
    for multiplier in multipliers:
        m = float(multiplier)
        if m <= 0.0 or not np.isfinite(m):
            raise ValueError("All bandwidth multipliers must be positive and finite.")
        scores[m] = loo_loglik_gaussian_kde(x, base * m)
    winner = max(scores, key=scores.get)
    return BandwidthSelection(
        base_scott=base,
        selected_multiplier=float(winner),
        selected_bandwidth=float(base * winner),
        scores=scores,
    )


def mode_normalizer(
    kde_centers: Sequence[float],
    bandwidth: float,
    k: float,
    a: float,
    b: float,
    gh_n: int = 40,
) -> float:
    """Compute Z = E_f0[exp(a cos(kX)+b sin(kX))] by Gauss-Hermite.

    For a Gaussian KDE, f0 is an equally weighted Gaussian mixture, so the
    expectation is evaluated component-by-component with deterministic
    Gauss-Hermite quadrature rather than Monte Carlo.
    """
    centers = np.asarray(kde_centers, dtype=float)
    if centers.size == 0:
        raise ValueError("KDE has no centers.")
    if bandwidth <= 0.0:
        raise ValueError("bandwidth must be positive")
    if gh_n < 8:
        raise ValueError("gh_n must be at least 8")
    nodes, weights = hermgauss(int(gh_n))
    x = centers[:, None] + math.sqrt(2.0) * bandwidth * nodes[None, :]
    r = a * np.cos(k * x) + b * np.sin(k * x)
    values = np.exp(np.clip(r, -50.0, 50.0))
    component_expectations = (values @ weights) / math.sqrt(math.pi)
    z = float(np.mean(component_expectations))
    if not np.isfinite(z) or z <= 0.0:
        raise RuntimeError("Non-finite residual-mode normalizer.")
    return z


def fixed_log_likelihood_ratio(
    x: Sequence[float], k: float, a: float, b: float, normalizer: float
) -> float:
    """Return 2 log[L(f1)/L(f0)] for a completely frozen residual model."""
    x = np.asarray(x, dtype=float)
    if normalizer <= 0.0 or not np.isfinite(normalizer):
        raise ValueError("normalizer must be positive and finite")
    residual = a * np.cos(k * x) + b * np.sin(k * x)
    return float(2.0 * np.sum(residual - math.log(normalizer)))


@dataclass(frozen=True)
class ModeAtK:
    k: float
    a: float
    b: float
    amplitude: float
    phase_atan2_b_a: float
    normalizer: float
    train_delta_2logl: float


def fit_mode_at_k(
    train_x: Sequence[float],
    kde_centers: Sequence[float],
    bandwidth: float,
    k: float,
    gh_n: int = 40,
) -> ModeAtK:
    train_x = np.asarray(train_x, dtype=float)

    def objective(theta: np.ndarray) -> float:
        a, b = float(theta[0]), float(theta[1])
        z = mode_normalizer(kde_centers, bandwidth, k, a, b, gh_n=gh_n)
        llr_half = np.sum(a * np.cos(k * train_x) + b * np.sin(k * train_x))
        llr_half -= train_x.size * math.log(z)
        return -float(llr_half)

    # Broad finite bounds protect the quadrature from pathological numerical
    # excursions. They are declared before holdout evaluation and are far above
    # the amplitudes seen in the V2 diagnostic.
    result = minimize(
        objective,
        np.zeros(2, dtype=float),
        method="L-BFGS-B",
        bounds=[(-3.0, 3.0), (-3.0, 3.0)],
    )
    if not result.success:
        raise RuntimeError(f"Mode fit failed at k={k}: {result.message}")
    a, b = map(float, result.x)
    z = mode_normalizer(kde_centers, bandwidth, k, a, b, gh_n=gh_n)
    delta = fixed_log_likelihood_ratio(train_x, k, a, b, z)
    return ModeAtK(
        k=float(k),
        a=a,
        b=b,
        amplitude=float(math.hypot(a, b)),
        phase_atan2_b_a=float(math.atan2(b, a)),
        normalizer=z,
        train_delta_2logl=delta,
    )


def scan_training_mode(
    train_x: Sequence[float],
    kde_centers: Sequence[float],
    bandwidth: float,
    k_grid: Sequence[float],
    gh_n: int = 40,
) -> tuple[ModeAtK, list[ModeAtK]]:
    results = [
        fit_mode_at_k(train_x, kde_centers, bandwidth, float(k), gh_n=gh_n)
        for k in k_grid
    ]
    best = max(results, key=lambda row: row.train_delta_2logl)
    return best, results


def sample_gaussian_kde(
    rng: np.random.Generator,
    kde_centers: Sequence[float],
    bandwidth: float,
    n: int,
) -> np.ndarray:
    centers = np.asarray(kde_centers, dtype=float)
    if n < 0:
        raise ValueError("n must be non-negative")
    index = rng.integers(0, centers.size, size=n)
    return centers[index] + bandwidth * rng.standard_normal(n)


def fixed_model_null(
    *,
    kde_centers: Sequence[float],
    bandwidth: float,
    n_events: int,
    k: float,
    a: float,
    b: float,
    normalizer: float,
    null_n: int,
    seed: int,
) -> np.ndarray:
    """Simulate the fixed holdout statistic under the frozen KDE baseline.

    Nothing is refit or rescanned in a null replicate.
    """
    if null_n <= 0:
        raise ValueError("null_n must be positive")
    rng = np.random.default_rng(seed)
    out = np.empty(null_n, dtype=float)
    for i in range(null_n):
        draw = sample_gaussian_kde(rng, kde_centers, bandwidth, n_events)
        out[i] = fixed_log_likelihood_ratio(draw, k, a, b, normalizer)
    return out
