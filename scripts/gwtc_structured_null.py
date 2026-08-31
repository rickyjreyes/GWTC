#!/usr/bin/env python3
"""Structured non-periodic nulls for the GWTC V4 robustness challenge.

The V4 null is deliberately more flexible than the featureless KDE null used
for the V3 fixed-model calibration.  In source-frame chirp mass z it uses a
continuous broken power law plus a truncated Gaussian peak, multiplied by a
monotone detection-efficiency weight,

    q_obs(z) propto q_intrinsic(z) * (z / z_max)**gamma.

The model contains no periodic term.  Parameters are fit on training data only.
The selection exponent gamma is not inferred; V4 treats several predeclared
values as sensitivity scenarios.  This is a phenomenological astrophysical
challenge, not a replacement for a full LVK population/selection analysis.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Sequence

import numpy as np
from scipy.optimize import minimize


@dataclass(frozen=True)
class StructuredNullParams:
    z_min: float
    z_max: float
    alpha_low: float
    alpha_high: float
    break_z: float
    peak_fraction: float
    peak_mu: float
    peak_sigma: float
    selection_gamma: float
    train_neg_loglik: float

    def to_dict(self) -> dict[str, float]:
        return {k: float(v) for k, v in asdict(self).items()}

    @classmethod
    def from_dict(cls, payload: dict[str, float]) -> "StructuredNullParams":
        return cls(**{k: float(v) for k, v in payload.items()})


def _validate_support(z_min: float, z_max: float) -> None:
    if not np.isfinite(z_min) or not np.isfinite(z_max):
        raise ValueError("support must be finite")
    if z_min <= 0.0 or z_max <= z_min:
        raise ValueError("require 0 < z_min < z_max")


def _grid(z_min: float, z_max: float, n: int) -> np.ndarray:
    _validate_support(z_min, z_max)
    if n < 256:
        raise ValueError("grid_n must be at least 256")
    # Chirp masses span multiplicative scales, so a geometric integration grid
    # resolves the low-mass region better than a uniform linear grid.
    return np.geomspace(z_min, z_max, int(n))


def broken_powerlaw_unnormalized(
    z: Sequence[float], alpha_low: float, alpha_high: float, break_z: float
) -> np.ndarray:
    """Continuous broken power law, before normalization."""
    x = np.asarray(z, dtype=float)
    if break_z <= 0.0 or not np.isfinite(break_z):
        raise ValueError("break_z must be positive and finite")
    ratio = np.clip(x / break_z, 1e-300, None)
    return np.where(x <= break_z, ratio ** alpha_low, ratio ** alpha_high)


def truncated_gaussian_unnormalized(
    z: Sequence[float], mu: float, sigma: float
) -> np.ndarray:
    x = np.asarray(z, dtype=float)
    if sigma <= 0.0 or not np.isfinite(sigma):
        raise ValueError("sigma must be positive and finite")
    u = (x - mu) / sigma
    return np.exp(-0.5 * np.clip(u * u, 0.0, 1400.0))


def selection_weight(
    z: Sequence[float], z_max: float, gamma: float
) -> np.ndarray:
    """Relative monotone efficiency used for the V4 sensitivity scenarios.

    For gamma=0 there is no selection weighting. gamma=2.5 corresponds to
    the familiar low-redshift inspiral-volume scaling M_chirp^(5/2), used only
    as an approximate stress test rather than an exact detector model.
    """
    x = np.asarray(z, dtype=float)
    if gamma < 0.0 or not np.isfinite(gamma):
        raise ValueError("selection_gamma must be finite and non-negative")
    if z_max <= 0.0:
        raise ValueError("z_max must be positive")
    return np.clip(x / z_max, 1e-300, 1.0) ** gamma


def observed_density_on_grid(
    grid: Sequence[float], params: StructuredNullParams
) -> np.ndarray:
    """Normalized selected density on an ordered support grid."""
    x = np.asarray(grid, dtype=float)
    if x.ndim != 1 or x.size < 2 or np.any(np.diff(x) <= 0.0):
        raise ValueError("grid must be one-dimensional and strictly increasing")
    if x[0] < params.z_min * (1.0 - 1e-12) or x[-1] > params.z_max * (1.0 + 1e-12):
        raise ValueError("grid lies outside parameter support")

    power = broken_powerlaw_unnormalized(
        x, params.alpha_low, params.alpha_high, params.break_z
    )
    peak = truncated_gaussian_unnormalized(x, params.peak_mu, params.peak_sigma)

    power_norm = float(np.trapezoid(power, x))
    peak_norm = float(np.trapezoid(peak, x))
    if power_norm <= 0.0 or peak_norm <= 0.0:
        raise RuntimeError("failed to normalize intrinsic components")
    power = power / power_norm
    peak = peak / peak_norm

    frac = float(params.peak_fraction)
    if not 0.0 <= frac <= 1.0:
        raise ValueError("peak_fraction must lie in [0,1]")
    intrinsic = (1.0 - frac) * power + frac * peak
    selected = intrinsic * selection_weight(x, params.z_max, params.selection_gamma)
    norm = float(np.trapezoid(selected, x))
    if not np.isfinite(norm) or norm <= 0.0:
        raise RuntimeError("failed to normalize selected structured null")
    return selected / norm


def observed_density_at(
    z: Sequence[float], params: StructuredNullParams, grid_n: int = 4096
) -> np.ndarray:
    """Evaluate the normalized selected density by interpolation."""
    x = np.asarray(z, dtype=float)
    grid = _grid(params.z_min, params.z_max, grid_n)
    density = observed_density_on_grid(grid, params)
    out = np.interp(x, grid, density, left=0.0, right=0.0)
    return np.asarray(out, dtype=float)


def _params_from_theta(
    theta: Sequence[float], z_min: float, z_max: float, gamma: float, nll: float
) -> StructuredNullParams:
    a1, a2, break_z, frac, mu, sigma = map(float, theta)
    return StructuredNullParams(
        z_min=float(z_min),
        z_max=float(z_max),
        alpha_low=a1,
        alpha_high=a2,
        break_z=break_z,
        peak_fraction=frac,
        peak_mu=mu,
        peak_sigma=sigma,
        selection_gamma=float(gamma),
        train_neg_loglik=float(nll),
    )


def fit_structured_null(
    z: Sequence[float], selection_gamma: float, grid_n: int = 2048
) -> StructuredNullParams:
    """Fit a non-periodic selected population to TRAINING masses only.

    A deterministic multi-start bounded optimization is used because the
    broken-power-law + peak likelihood is mildly non-convex.
    """
    data = np.asarray(z, dtype=float)
    data = data[np.isfinite(data) & (data > 0.0)]
    if data.size < 30:
        raise ValueError("Need at least 30 positive finite training values.")
    if selection_gamma < 0.0 or not np.isfinite(selection_gamma):
        raise ValueError("selection_gamma must be finite and non-negative")

    lo_data = float(np.min(data))
    hi_data = float(np.max(data))
    eps = 1e-8 * (hi_data - lo_data)
    z_min = max(np.nextafter(0.0, 1.0), lo_data - eps)
    z_max = hi_data + eps
    span = z_max - z_min
    grid = _grid(z_min, z_max, grid_n)

    break_lo = z_min + 0.10 * span
    break_hi = z_max - 0.10 * span
    sigma_lo = max(1e-6, 0.015 * span)
    sigma_hi = 0.50 * span

    bounds = [
        (-8.0, 4.0),
        (-10.0, 4.0),
        (break_lo, break_hi),
        (0.0, 0.75),
        (z_min, z_max),
        (sigma_lo, sigma_hi),
    ]

    q25, q50, q75 = np.quantile(data, [0.25, 0.50, 0.75])
    scale = max(float(np.std(data, ddof=1)), sigma_lo)
    starts = [
        [-1.5, -3.0, q50, 0.20, q75, min(max(0.20 * scale, sigma_lo), sigma_hi)],
        [-2.0, -4.0, q75, 0.10, q50, min(max(0.35 * scale, sigma_lo), sigma_hi)],
        [-0.5, -2.5, q25, 0.30, q75, min(max(0.15 * span, sigma_lo), sigma_hi)],
        [-3.0, -1.0, q50, 0.05, q25, min(max(0.25 * span, sigma_lo), sigma_hi)],
    ]

    def objective(theta: np.ndarray) -> float:
        provisional = _params_from_theta(theta, z_min, z_max, selection_gamma, 0.0)
        try:
            density_grid = observed_density_on_grid(grid, provisional)
        except (ValueError, RuntimeError, FloatingPointError):
            return 1e100
        dens = np.interp(data, grid, density_grid, left=0.0, right=0.0)
        dens = np.clip(dens, 1e-300, None)
        value = -float(np.sum(np.log(dens)))
        return value if np.isfinite(value) else 1e100

    winners: list[tuple[float, np.ndarray]] = []
    for start in starts:
        clipped = np.asarray(
            [min(max(v, low), high) for v, (low, high) in zip(start, bounds)],
            dtype=float,
        )
        result = minimize(
            objective,
            clipped,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 2000, "ftol": 1e-11},
        )
        if result.success and np.isfinite(result.fun):
            winners.append((float(result.fun), np.asarray(result.x, dtype=float)))

    if not winners:
        raise RuntimeError("All structured-null optimization starts failed.")
    best_nll, best_theta = min(winners, key=lambda item: item[0])
    return _params_from_theta(
        best_theta, z_min, z_max, selection_gamma, best_nll
    )


def sample_structured_null(
    rng: np.random.Generator,
    params: StructuredNullParams,
    n: int,
    grid_n: int = 4096,
) -> np.ndarray:
    """Draw from the fitted selected density by deterministic-grid inversion."""
    if n < 0:
        raise ValueError("n must be non-negative")
    grid = _grid(params.z_min, params.z_max, grid_n)
    density = observed_density_on_grid(grid, params)
    dx = np.diff(grid)
    increments = 0.5 * (density[:-1] + density[1:]) * dx
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    if cdf[-1] <= 0.0 or not np.isfinite(cdf[-1]):
        raise RuntimeError("invalid structured-null CDF")
    cdf /= cdf[-1]
    u = rng.random(int(n))
    return np.interp(u, cdf, grid)


def structured_fixed_statistic_null(
    *,
    params: StructuredNullParams,
    n_events: int,
    k: float,
    a: float,
    b: float,
    normalizer: float,
    null_n: int,
    seed: int,
    grid_n: int = 4096,
) -> np.ndarray:
    """Monte Carlo distribution of the completely fixed residual statistic.

    No population parameters, residual coefficients, frequency, or bandwidth
    are fit inside a replicate.
    """
    if null_n <= 0:
        raise ValueError("null_n must be positive")
    if n_events <= 0:
        raise ValueError("n_events must be positive")

    from gwtc_unbinned_kde import fixed_log_likelihood_ratio

    rng = np.random.default_rng(int(seed))
    out = np.empty(int(null_n), dtype=float)
    for i in range(int(null_n)):
        z = sample_structured_null(rng, params, int(n_events), grid_n=grid_n)
        out[i] = fixed_log_likelihood_ratio(
            np.log(z), float(k), float(a), float(b), float(normalizer)
        )
    return out
