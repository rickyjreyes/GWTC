from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_baseline_v2_cv.py"
spec = importlib.util.spec_from_file_location("baseline_v2", SCRIPT)
assert spec is not None and spec.loader is not None
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_poisson_deviance_is_zero_on_exact_mean() -> None:
    y = np.array([0.0, 1.0, 3.0, 8.0])
    assert abs(mod.poisson_deviance(y, y)) < 1e-10


def test_polynomial_poisson_fit_recovers_smooth_counts() -> None:
    x = np.linspace(-1.0, 1.0, 40)
    center = float(np.mean(x))
    scale = float(np.std(x))
    X = mod.design_poly(x, degree=2, center=center, scale=scale)
    beta_true = np.array([2.0, 0.25, -0.35])
    mu = np.exp(X @ beta_true)
    # Fit the exact expected counts rather than a random draw so this is deterministic.
    beta = mod.fit_poisson_glm(X, mu)
    mu_fit = np.exp(X @ beta)
    assert mod.poisson_deviance(mu, mu_fit) < 1e-6


def test_folds_cover_each_bin_once() -> None:
    folds = mod.make_folds(17, 5)
    merged = np.concatenate(folds)
    assert sorted(merged.tolist()) == list(range(17))
