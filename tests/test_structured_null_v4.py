from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from freeze_gwtc_v4_structured_null import K_BAND, K_PRED  # noqa: E402
from gwtc_structured_null import (  # noqa: E402
    StructuredNullParams,
    fit_structured_null,
    observed_density_on_grid,
    sample_structured_null,
    selection_weight,
    structured_fixed_statistic_null,
)
from gwtc_unbinned_kde import mode_normalizer  # noqa: E402


def synthetic_params(gamma: float = 1.5) -> StructuredNullParams:
    return StructuredNullParams(
        z_min=5.0,
        z_max=60.0,
        alpha_low=-1.2,
        alpha_high=-3.0,
        break_z=24.0,
        peak_fraction=0.22,
        peak_mu=32.0,
        peak_sigma=4.5,
        selection_gamma=gamma,
        train_neg_loglik=0.0,
    )


def test_published_frequency_is_hard_fixed() -> None:
    assert K_PRED == 9.7
    assert K_BAND == (9.5, 10.0)


def test_frozen_v4_is_phase_specific_and_null_is_nonperiodic() -> None:
    """Protect the scientific separation built into the committed V4 freeze.

    The signal statistic must keep a nonzero, finite training-frozen
    amplitude/phase at exactly k=9.7. The structured population scenarios are
    allowed broad non-periodic mass structure and selection weighting, but no
    periodic/log-periodic parameter may be inserted into their parameter sets.
    """
    path = ROOT / "tables" / "gwtc_v4_frozen_structured_null.json"
    payload = json.loads(path.read_text(encoding="utf-8"))

    assert payload["schema"] == "GWTC_V4_STRUCTURED_NULL_FIXED_PREDICTION_V1"
    assert payload["prediction_status"] == "EXTERNALLY_PUBLISHED_AND_FROZEN_BEFORE_V4"
    assert payload["holdout_evaluated_by_this_script"] is False
    assert payload["k_pred"] == 9.7

    mode = payload["fixed_mode"]
    assert mode["k"] == 9.7
    assert np.isfinite(mode["amplitude"])
    assert mode["amplitude"] > 0.0
    assert np.isfinite(mode["phase_atan2_b_a"])
    assert np.isfinite(mode["a_cos"])
    assert np.isfinite(mode["b_sin"])
    np.testing.assert_allclose(
        mode["amplitude"],
        np.hypot(mode["a_cos"], mode["b_sin"]),
        rtol=1e-12,
        atol=1e-12,
    )

    scenarios = payload["structured_null_scenarios"]
    assert [row["selection_gamma"] for row in scenarios] == [0.0, 1.5, 2.5]

    forbidden_fragments = ("period", "frequency", "phase", "cos", "sin", "k_mode")
    for scenario in scenarios:
        assert scenario["population_family"] == (
            "continuous_broken_powerlaw_plus_truncated_gaussian_peak"
        )
        params = scenario["params"]
        lowered_keys = [str(key).lower() for key in params]
        assert not any(
            fragment in key
            for key in lowered_keys
            for fragment in forbidden_fragments
        )
        assert set(params) == {
            "z_min",
            "z_max",
            "alpha_low",
            "alpha_high",
            "break_z",
            "peak_fraction",
            "peak_mu",
            "peak_sigma",
            "selection_gamma",
            "train_neg_loglik",
        }


def test_selection_weight_gamma_zero_is_unity() -> None:
    z = np.array([5.0, 20.0, 60.0])
    w = selection_weight(z, 60.0, 0.0)
    np.testing.assert_allclose(w, np.ones_like(z))


def test_structured_density_integrates_to_one() -> None:
    params = synthetic_params()
    grid = np.geomspace(params.z_min, params.z_max, 5000)
    density = observed_density_on_grid(grid, params)
    integral = np.trapezoid(density, grid)
    assert np.all(np.isfinite(density))
    assert np.all(density >= 0.0)
    assert abs(integral - 1.0) < 2e-5


def test_structured_sampling_is_deterministic_and_inside_support() -> None:
    params = synthetic_params()
    a = sample_structured_null(np.random.default_rng(1234), params, 200, grid_n=1024)
    b = sample_structured_null(np.random.default_rng(1234), params, 200, grid_n=1024)
    np.testing.assert_allclose(a, b, rtol=0.0, atol=0.0)
    assert np.min(a) >= params.z_min
    assert np.max(a) <= params.z_max


def test_structured_fit_is_training_only_deterministic() -> None:
    source = synthetic_params(gamma=1.5)
    train = sample_structured_null(np.random.default_rng(44), source, 90, grid_n=1024)
    first = fit_structured_null(train, 1.5, grid_n=512)
    second = fit_structured_null(train, 1.5, grid_n=512)
    for key, value in first.to_dict().items():
        assert np.isfinite(value)
        assert abs(value - second.to_dict()[key]) < 1e-9
    assert first.z_min <= float(np.min(train))
    assert first.z_max >= float(np.max(train))
    assert 0.0 <= first.peak_fraction <= 0.75


def test_structured_fixed_statistic_null_is_deterministic() -> None:
    params = synthetic_params(gamma=2.5)
    centers = np.linspace(np.log(params.z_min), np.log(params.z_max), 35)
    bandwidth = 0.18
    k = 9.7
    a = 0.08
    b = 0.22
    z_norm = mode_normalizer(centers, bandwidth, k, a, b, gh_n=20)

    first = structured_fixed_statistic_null(
        params=params,
        n_events=25,
        k=k,
        a=a,
        b=b,
        normalizer=z_norm,
        null_n=30,
        seed=999,
        grid_n=512,
    )
    second = structured_fixed_statistic_null(
        params=params,
        n_events=25,
        k=k,
        a=a,
        b=b,
        normalizer=z_norm,
        null_n=30,
        seed=999,
        grid_n=512,
    )
    np.testing.assert_allclose(first, second, rtol=0.0, atol=0.0)
    assert first.shape == (30,)
    assert np.all(np.isfinite(first))
