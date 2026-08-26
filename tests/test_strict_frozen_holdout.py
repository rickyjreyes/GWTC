from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import evaluate_frozen_gwtc_holdout as E  # noqa: E402
import fit_frozen_gwtc_mode as F  # noqa: E402
import make_holdout_manifest as H  # noqa: E402


def test_manifest_uses_only_canonical_primary_entries() -> None:
    df = pd.DataFrame(
        {
            "commonName": ["A", "A", "B", "C"],
            "catalog_version": ["GWTC-4.1", "GWTC-5.0", "GWTC-5.0", "GWTC-4.1"],
            "p_astro": [0.99, 0.99, 0.95, 0.9],
            "is_primary_entry": [True, False, True, True],
        }
    )
    manifest = H.build_manifest(df, holdout_prefix="GWTC-5", p_astro_min=0.5)
    got = manifest.set_index("event_name")
    assert len(got) == 3
    assert got.loc["A", "split"] == "train"
    assert got.loc["B", "split"] == "holdout"
    assert got.loc["C", "split"] == "train"


def test_predictive_delta_prefers_matching_frozen_model() -> None:
    counts = np.array([5, 15, 80])
    p0 = np.array([1 / 3, 1 / 3, 1 / 3])
    p1 = np.array([0.05, 0.15, 0.80])
    assert E.predictive_delta(counts, p0, p1) > 0
    assert E.predictive_delta(counts, p1, p0) < 0


def test_fixed_model_null_is_deterministic_and_detects_strong_prediction() -> None:
    counts = np.array([5, 15, 80])
    p0 = np.array([1 / 3, 1 / 3, 1 / 3])
    p1 = np.array([0.05, 0.15, 0.80])
    p_a, null_a = E.fixed_model_null_p(counts, p0, p1, null_n=500, seed=7)
    p_b, null_b = E.fixed_model_null_p(counts, p0, p1, null_n=500, seed=7)
    assert p_a == p_b
    assert np.array_equal(null_a, null_b)
    assert p_a < 0.05


def test_training_fit_freezes_probabilities_and_respects_nyquist() -> None:
    rng = np.random.default_rng(42)
    ell = rng.normal(3.0, 0.35, size=800)
    values = np.exp(ell)
    model = F.fit_training_mode(
        values,
        baseline_degree=4,
        bins=40,
        k_min=0.5,
        k_max=1e6,
        n_k=60,
    )
    p0 = np.asarray(model["baseline_probabilities"])
    p1 = np.asarray(model["residual_probabilities"])
    assert np.isclose(p0.sum(), 1.0)
    assert np.isclose(p1.sum(), 1.0)
    assert model["k_grid_effective_max"] <= model["nyquist_k"] + 1e-12
    assert model["k_star"] <= model["nyquist_k"] + 1e-12
    assert model["holdout_evaluated"] is False
