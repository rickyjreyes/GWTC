#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
wct_gwtc_lib.py
===============

Core library for the WCT-motivated log-domain residual diagnostic on the
LIGO-Virgo-KAGRA Gravitational-Wave Transient Catalog (GWTC).

This is a DIAGNOSTIC HARNESS. It does not prove WCT. It does not replace LVK
population inference. It only tests whether selected catalog-scale variables
contain a stable log-domain residual mode that survives null testing and
stress checks.

Core model
----------
For a positive, physically scale-like catalog variable z > 0, define the
log coordinate:

    ell_i = ln(z_i)

Bin ell into B bins to form a 1-D event-density (count) field y_j over bin
centers ell_j. Fit a smooth Poisson baseline mu0(ell) (low-order polynomial
in log-rate). Then test a single log-periodic residual mode:

    log mu(ell; k, a, b, c) = log mu0(ell) + c + a cos(k ell) + b sin(k ell)

For each k in a declared grid, fit (a, b, c) by Poisson regression and record
the deviance improvement:

    DeltaD(k) = D[y || mu0] - D[y || mu(k)]

Primary statistic:

    T_obs = max_k DeltaD(k)
    k_star = argmax_k DeltaD(k)
    n_star = k_star * Delta_ell_A / (2 * pi)

where Delta_ell_A is the retained active log-domain support (ell_max - ell_min).

Null model
----------
Parametric Poisson bootstrap from the fitted baseline mu0:
    1. draw y_null ~ Poisson(mu0)
    2. refit the same baseline
    3. scan the same k-grid
    4. record T_null = max_k DeltaD_null(k)

Global p-value:

    p_global = (1 + #{T_null >= T_obs}) / (1 + N_null)

PASS is NEVER declared from a local chi^2 p-value alone.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Admissible scale-like variables -> column in the clean events table.
# These are physically scale-like (positive, ratio-meaningful) variables.
SCALE_VARIABLES = {
    "M_chirp": "M_chirp",
    "M_total": "M_total",
    "D_L": "D_L",
    "M_final": "M_final",
    "redshift": "redshift",
}

# Bounded / non-scale variables: DIAGNOSTIC ONLY. A primary PASS may never be
# declared from these. chi_eff is bounded in [-1, 1]; we expose |chi_eff| as a
# deliberate negative-control style coordinate, not a scale test.
BOUNDED_VARIABLES = {
    "abs_chi_eff": "chi_eff",
}

# ---------------------------------------------------------------------------
# Canonical CSV output schema. Every result CSV in this harness uses these
# columns (extra columns may be appended after).
# ---------------------------------------------------------------------------
CSV_COLUMNS: List[str] = [
    "catalog_version",
    "variable",
    "subset",
    "bin_count",
    "k_best",
    "n_star",
    "DeltaD_star",
    "local_p",
    "global_p",
    "null_n",
    "verdict_label",
    "notes",
]

# Allowed verdict labels for the per-result column.
VERDICT_LABELS = ("PASS", "PARTIAL", "FAIL", "INCOMPLETE")

# Reliability classes.
RELIABILITY_CLASSES = ("I", "II", "III", "IV")


# ---------------------------------------------------------------------------
# Poisson deviance and GLM fitting (IRLS)
# ---------------------------------------------------------------------------
def poisson_deviance(y: np.ndarray, mu: np.ndarray) -> float:
    """Poisson deviance D = 2 * sum[ y*log(y/mu) - (y - mu) ].

    Uses the convention 0*log(0/mu) = 0.
    """
    y = np.asarray(y, dtype=float)
    mu = np.asarray(mu, dtype=float)
    mu = np.clip(mu, 1e-12, None)
    term = np.zeros_like(y)
    nz = y > 0
    term[nz] = y[nz] * np.log(y[nz] / mu[nz])
    dev = 2.0 * np.sum(term - (y - mu))
    return float(dev)


def fit_poisson_glm(
    y: np.ndarray,
    X: np.ndarray,
    offset: Optional[np.ndarray] = None,
    max_iter: int = 100,
    tol: float = 1e-9,
) -> Tuple[np.ndarray, np.ndarray]:
    """Fit a Poisson GLM with log link by iteratively reweighted least squares.

    Model: log(mu) = offset + X @ beta

    Returns (beta, mu). Robust to rank-deficiency via a small ridge.
    """
    y = np.asarray(y, dtype=float)
    X = np.asarray(X, dtype=float)
    n, p = X.shape
    if offset is None:
        offset = np.zeros(n)
    offset = np.asarray(offset, dtype=float)

    beta = np.zeros(p)
    # Sensible start for intercept-like first column.
    eta = offset + X @ beta
    mu = np.exp(np.clip(eta, -30, 30))

    ridge = 1e-8
    for _ in range(max_iter):
        mu = np.clip(mu, 1e-12, None)
        # IRLS working response and weights for Poisson/log link:
        # W = mu, z = eta - offset + (y - mu) / mu
        W = mu
        z = (eta - offset) + (y - mu) / mu
        WX = X * W[:, None]
        A = X.T @ WX + ridge * np.eye(p)
        bvec = X.T @ (W * z)
        try:
            beta_new = np.linalg.solve(A, bvec)
        except np.linalg.LinAlgError:
            beta_new = np.linalg.lstsq(A, bvec, rcond=None)[0]
        eta = offset + X @ beta_new
        eta = np.clip(eta, -30, 30)
        mu_new = np.exp(eta)
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            mu = mu_new
            break
        beta = beta_new
        mu = mu_new
    return beta, mu


def polynomial_design(ell: np.ndarray, degree: int) -> np.ndarray:
    """Design matrix [1, t, t^2, ..., t^degree] on centered/scaled ell.

    ell is centered and scaled to unit std for numerical stability.
    """
    ell = np.asarray(ell, dtype=float)
    c = np.mean(ell)
    s = np.std(ell)
    if s == 0:
        s = 1.0
    t = (ell - c) / s
    cols = [np.ones_like(t)]
    for d in range(1, degree + 1):
        cols.append(t ** d)
    return np.column_stack(cols)


# ---------------------------------------------------------------------------
# Histogram / event-density construction
# ---------------------------------------------------------------------------
@dataclass
class DensityField:
    ell: np.ndarray              # event log values used
    edges: np.ndarray            # bin edges
    centers: np.ndarray          # bin centers
    counts: np.ndarray           # y_j
    ell_min: float
    ell_max: float
    delta_ell_active: float      # retained active support = ell_max - ell_min
    bin_count: int


def build_density_field(
    z: Sequence[float],
    bin_count: int,
    ell_min: Optional[float] = None,
    ell_max: Optional[float] = None,
) -> DensityField:
    """Build a binned 1-D event-density field of ell = ln(z).

    Only strictly positive, finite z values are retained.
    """
    z = np.asarray(z, dtype=float)
    z = z[np.isfinite(z) & (z > 0)]
    if z.size == 0:
        raise ValueError("No positive finite values to build density field.")
    ell = np.log(z)
    lo = float(np.min(ell)) if ell_min is None else float(ell_min)
    hi = float(np.max(ell)) if ell_max is None else float(ell_max)
    if hi <= lo:
        hi = lo + 1e-6
    edges = np.linspace(lo, hi, bin_count + 1)
    counts, _ = np.histogram(ell, bins=edges)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return DensityField(
        ell=ell,
        edges=edges,
        centers=centers,
        counts=counts.astype(float),
        ell_min=lo,
        ell_max=hi,
        delta_ell_active=hi - lo,
        bin_count=bin_count,
    )


# ---------------------------------------------------------------------------
# Baseline + log-periodic scan
# ---------------------------------------------------------------------------
def fit_baseline(field: DensityField, degree: int = 3) -> np.ndarray:
    """Fit smooth Poisson baseline mu0(ell) via polynomial log-rate GLM."""
    X = polynomial_design(field.centers, degree)
    _, mu0 = fit_poisson_glm(field.counts, X)
    return mu0


def delta_d_at_k(
    counts: np.ndarray,
    centers: np.ndarray,
    mu0: np.ndarray,
    k: float,
) -> float:
    """Deviance improvement for one log-periodic mode at wavenumber k.

    log mu = log mu0 + c + a cos(k ell) + b sin(k ell)
    """
    offset = np.log(np.clip(mu0, 1e-12, None))
    X = np.column_stack([
        np.ones_like(centers),
        np.cos(k * centers),
        np.sin(k * centers),
    ])
    _, mu = fit_poisson_glm(counts, X, offset=offset)
    d0 = poisson_deviance(counts, mu0)
    dk = poisson_deviance(counts, mu)
    return d0 - dk


@dataclass
class ScanResult:
    k_grid: np.ndarray
    delta_d: np.ndarray
    k_star: float
    delta_d_star: float
    n_star: float
    delta_ell_active: float
    bin_count: int
    local_p: Optional[float] = None


def scan_k(
    field: DensityField,
    k_grid: Sequence[float],
    baseline_degree: int = 3,
    mu0: Optional[np.ndarray] = None,
) -> ScanResult:
    """Scan the k-grid and return the best log-periodic mode.

    n_star = k_star * Delta_ell_A / (2*pi).
    A local (asymptotic) p-value is computed from the chi^2_2 tail of the
    deviance improvement at k_star. This local p-value is DIAGNOSTIC ONLY and
    must never be used on its own to declare PASS.
    """
    k_grid = np.asarray(list(k_grid), dtype=float)
    if mu0 is None:
        mu0 = fit_baseline(field, degree=baseline_degree)
    dd = np.array([
        delta_d_at_k(field.counts, field.centers, mu0, k) for k in k_grid
    ])
    idx = int(np.argmax(dd))
    k_star = float(k_grid[idx])
    dd_star = float(dd[idx])
    n_star = k_star * field.delta_ell_active / (2.0 * math.pi)

    # Local p-value: chi^2 with 2 dof (a, b added beyond the constant c).
    # This is the *uncorrected, local* tail; it ignores the scan over k.
    from scipy.stats import chi2
    local_p = float(chi2.sf(max(dd_star, 0.0), df=2))

    return ScanResult(
        k_grid=k_grid,
        delta_d=dd,
        k_star=k_star,
        delta_d_star=dd_star,
        n_star=n_star,
        delta_ell_active=field.delta_ell_active,
        bin_count=field.bin_count,
        local_p=local_p,
    )


# ---------------------------------------------------------------------------
# Parametric Poisson bootstrap null
# ---------------------------------------------------------------------------
@dataclass
class NullResult:
    t_obs: float
    t_null: np.ndarray
    global_p: float
    null_n: int


def run_poisson_null(
    field: DensityField,
    k_grid: Sequence[float],
    t_obs: float,
    null_n: int,
    seed: int = 12345,
    baseline_degree: int = 3,
    refit_baseline: bool = True,
) -> NullResult:
    """Parametric Poisson bootstrap from the fitted baseline mu0.

    For each replicate: draw y ~ Poisson(mu0), (optionally) refit the same
    baseline, scan the same k-grid, record T_null = max_k DeltaD.

    p_global = (1 + #{T_null >= T_obs}) / (1 + N_null)
    """
    rng = np.random.default_rng(seed)
    mu0 = fit_baseline(field, degree=baseline_degree)
    k_grid = np.asarray(list(k_grid), dtype=float)
    t_null = np.empty(null_n)
    for j in range(null_n):
        y = rng.poisson(mu0).astype(float)
        if refit_baseline:
            Xb = polynomial_design(field.centers, baseline_degree)
            _, mu0_j = fit_poisson_glm(y, Xb)
        else:
            mu0_j = mu0
        dd = np.array([
            delta_d_at_k(y, field.centers, mu0_j, k) for k in k_grid
        ])
        t_null[j] = float(np.max(dd))
    n_ge = int(np.sum(t_null >= t_obs))
    p_global = (1 + n_ge) / (1 + null_n)
    return NullResult(t_obs=t_obs, t_null=t_null, global_p=p_global, null_n=null_n)


def global_p_value(t_null: Sequence[float], t_obs: float) -> float:
    """p_global = (1 + #{T_null >= T_obs}) / (1 + N_null). Pure formula."""
    t_null = np.asarray(t_null, dtype=float)
    n = t_null.size
    n_ge = int(np.sum(t_null >= t_obs))
    return (1 + n_ge) / (1 + n)


def compute_n_star(k_star: float, delta_ell_active: float) -> float:
    """n_star = k_star * Delta_ell_A / (2*pi). Pure formula."""
    return k_star * delta_ell_active / (2.0 * math.pi)


# ---------------------------------------------------------------------------
# Verdict logic (per-result heuristic; the master verdict is assembled by
# make_gwtc_verdict.py using the full stress-test matrix).
# ---------------------------------------------------------------------------
def per_result_verdict(
    global_p: Optional[float],
    null_n: int,
    have_data: bool = True,
    min_null_n: int = 1000,
) -> str:
    """Heuristic single-result verdict label.

    This is intentionally conservative. A real PASS requires the full
    stress-test matrix (bin / variable / threshold / run / controls) which is
    adjudicated by make_gwtc_verdict.py, not here.
    """
    if not have_data or global_p is None:
        return "INCOMPLETE"
    if null_n < min_null_n:
        # Not enough nulls to support a PASS-grade claim under this harness.
        return "PARTIAL" if global_p < 0.05 else "FAIL"
    if global_p < 0.05:
        return "PASS"
    return "FAIL"


def make_k_grid(k_min: float, k_max: float, n_k: int) -> np.ndarray:
    """Uniform k grid. Declared explicitly so the search family is auditable."""
    return np.linspace(k_min, k_max, n_k)


def result_row(
    catalog_version: str,
    variable: str,
    subset: str,
    bin_count: int,
    k_best: float,
    n_star: float,
    delta_d_star: float,
    local_p: Optional[float],
    global_p: Optional[float],
    null_n: int,
    verdict_label: str,
    notes: str = "",
) -> Dict[str, object]:
    """Build a canonical result-row dict matching CSV_COLUMNS."""
    return {
        "catalog_version": catalog_version,
        "variable": variable,
        "subset": subset,
        "bin_count": bin_count,
        "k_best": k_best,
        "n_star": n_star,
        "DeltaD_star": delta_d_star,
        "local_p": local_p if local_p is not None else "",
        "global_p": global_p if global_p is not None else "",
        "null_n": null_n,
        "verdict_label": verdict_label,
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# Data loading and subset selection
# ---------------------------------------------------------------------------
def load_clean_table(path: str = os.path.join("tables", "gwtc_events_clean.csv")):
    """Load the clean per-event table built by build_gwtc_table.py."""
    import pandas as pd
    return pd.read_csv(path)


def select_subset(
    df,
    catalog: Optional[str] = None,
    cumulative: bool = False,
    source_class: Optional[str] = None,
    p_astro_min: Optional[float] = None,
    far_max: Optional[float] = None,
    observing_run: Optional[str] = None,
):
    """Apply selection cuts and return a filtered copy.

    catalog       : keep one catalog_version (e.g. 'GWTC-5.0').
    cumulative    : keep only deduplicated primary entries (one per event).
    source_class  : 'BBH' | 'NSBH' | 'BNS'.
    p_astro_min   : keep rows with p_astro >= threshold.
    far_max       : keep rows with far <= threshold (per year).
    observing_run : keep one run label (e.g. 'O4a').
    """
    out = df.copy()
    if cumulative:
        out = out[out["is_primary_entry"] == True]  # noqa: E712
    if catalog is not None:
        out = out[out["catalog_version"] == catalog]
    if source_class is not None:
        out = out[out["source_class"] == source_class]
    if p_astro_min is not None:
        out = out[out["p_astro"].fillna(-1) >= p_astro_min]
    if far_max is not None:
        out = out[out["far"].notna() & (out["far"] <= far_max)]
    if observing_run is not None:
        out = out[out["observing_run"] == observing_run]
    return out


def get_variable_values(df, variable: str):
    """Return finite, positive values for an admissible scale variable, or the
    bounded-coordinate transform for diagnostic variables.

    For bounded 'abs_chi_eff' we return |chi_eff|, which is NOT a scale
    variable and is used only as a coordinate-stress / negative-control test.
    """
    if variable in SCALE_VARIABLES:
        col = SCALE_VARIABLES[variable]
        vals = df[col].to_numpy(dtype=float)
        vals = vals[np.isfinite(vals) & (vals > 0)]
        return vals
    if variable in BOUNDED_VARIABLES:
        col = BOUNDED_VARIABLES[variable]
        vals = np.abs(df[col].to_numpy(dtype=float))
        vals = vals[np.isfinite(vals) & (vals > 0)]
        return vals
    raise KeyError(f"Unknown variable '{variable}'.")


# ---------------------------------------------------------------------------
# One-cell analysis: scan + null + verdict, returning a canonical result row.
# This is the shared code path used by every stress script so behaviour is
# identical across the harness.
# ---------------------------------------------------------------------------
def analyze_cell(
    values: Sequence[float],
    variable: str,
    catalog_version: str,
    subset: str,
    bin_count: int,
    k_grid: Sequence[float],
    null_n: int,
    seed: int = 12345,
    baseline_degree: int = 3,
    min_events: int = 20,
    min_null_n: int = 1000,
    notes: str = "",
) -> Dict[str, object]:
    """Run the full single-cell diagnostic and return a canonical result row.

    Returns INCOMPLETE if there are too few events to bin meaningfully.
    Bounded (non-scale) variables are flagged in notes and can never be a
    primary PASS at the master-verdict stage.
    """
    values = np.asarray(list(values), dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    n_events = int(values.size)
    bounded = variable in BOUNDED_VARIABLES
    note_parts = [p for p in [notes] if p]
    if bounded:
        note_parts.append("bounded_non_scale_diagnostic_only")

    if n_events < min_events:
        note_parts.append(f"n_events={n_events}<min_events={min_events}")
        return result_row(
            catalog_version, variable, subset, bin_count,
            k_best=float("nan"), n_star=float("nan"), delta_d_star=float("nan"),
            local_p=None, global_p=None, null_n=0,
            verdict_label="INCOMPLETE", notes=";".join(note_parts),
        )

    field = build_density_field(values, bin_count=bin_count)
    mu0 = fit_baseline(field, degree=baseline_degree)
    scan = scan_k(field, k_grid, baseline_degree=baseline_degree, mu0=mu0)

    global_p = None
    if null_n > 0:
        null = run_poisson_null(
            field, k_grid, scan.delta_d_star, null_n=null_n,
            seed=seed, baseline_degree=baseline_degree, refit_baseline=True,
        )
        global_p = null.global_p

    verdict = per_result_verdict(global_p, null_n, have_data=True, min_null_n=min_null_n)
    note_parts.append(f"n_events={n_events}")
    return result_row(
        catalog_version, variable, subset, bin_count,
        k_best=scan.k_star, n_star=scan.n_star, delta_d_star=scan.delta_d_star,
        local_p=scan.local_p, global_p=global_p, null_n=null_n,
        verdict_label=verdict, notes=";".join(note_parts),
    )
