from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from gwtc_unbinned_kde import (  # noqa: E402
    fixed_log_likelihood_ratio,
    fixed_model_null,
    fit_mode_at_k,
    mode_normalizer,
    select_bandwidth_loo,
)


def test_zero_mode_normalizer_is_one() -> None:
    centers = np.array([-1.0, -0.2, 0.4, 1.3])
    z = mode_normalizer(centers, 0.25, k=7.0, a=0.0, b=0.0, gh_n=24)
    assert abs(z - 1.0) < 1e-12


def test_bandwidth_selection_is_finite_and_declared() -> None:
    rng = np.random.default_rng(12)
    x = rng.normal(size=80)
    selection = select_bandwidth_loo(x, [0.7, 1.0, 1.4])
    assert selection.selected_multiplier in {0.7, 1.0, 1.4}
    assert selection.selected_bandwidth > 0.0
    assert all(np.isfinite(list(selection.scores.values())))


def test_mode_fit_cannot_be_worse_than_zero_tilt() -> None:
    rng = np.random.default_rng(3)
    x = np.concatenate([rng.normal(-0.7, 0.18, 45), rng.normal(0.7, 0.18, 45)])
    mode = fit_mode_at_k(x, x, bandwidth=0.25, k=4.0, gh_n=24)
    assert np.isfinite(mode.train_delta_2logl)
    assert mode.train_delta_2logl >= -1e-7
    assert mode.amplitude >= 0.0


def test_fixed_model_null_is_deterministic() -> None:
    centers = np.array([-0.8, -0.1, 0.5, 1.0])
    z = mode_normalizer(centers, 0.3, k=3.0, a=0.2, b=-0.1, gh_n=24)
    a = fixed_model_null(
        kde_centers=centers,
        bandwidth=0.3,
        n_events=12,
        k=3.0,
        a=0.2,
        b=-0.1,
        normalizer=z,
        null_n=50,
        seed=99,
    )
    b = fixed_model_null(
        kde_centers=centers,
        bandwidth=0.3,
        n_events=12,
        k=3.0,
        a=0.2,
        b=-0.1,
        normalizer=z,
        null_n=50,
        seed=99,
    )
    np.testing.assert_allclose(a, b, rtol=0.0, atol=0.0)


def test_fixed_llr_uses_no_baseline_refit() -> None:
    x = np.array([0.0, 0.1, 0.2])
    got = fixed_log_likelihood_ratio(x, k=2.0, a=0.3, b=0.0, normalizer=1.1)
    expected = 2.0 * np.sum(0.3 * np.cos(2.0 * x) - np.log(1.1))
    assert abs(got - expected) < 1e-12


def test_freeze_cli_reads_training_split_only(tmp_path: Path) -> None:
    rng = np.random.default_rng(44)
    train_n = 60
    holdout_n = 40
    train_names = [f"T{i:03d}" for i in range(train_n)]
    holdout_names = [f"H{i:03d}" for i in range(holdout_n)]

    # Holdout is deliberately far away. If it leaks into fitting, the frozen
    # training span and training count will expose it.
    train_z = np.exp(rng.normal(2.0, 0.15, train_n))
    holdout_z = np.exp(rng.normal(5.0, 0.15, holdout_n))
    table = pd.DataFrame(
        {
            "commonName": train_names + holdout_names,
            "is_primary_entry": [True] * (train_n + holdout_n),
            "M_chirp": np.concatenate([train_z, holdout_z]),
            "p_astro": [0.99] * (train_n + holdout_n),
        }
    )
    manifest = pd.DataFrame(
        {
            "event_name": train_names + holdout_names,
            "split": ["train"] * train_n + ["holdout"] * holdout_n,
        }
    )

    table_path = tmp_path / "events.csv"
    manifest_path = tmp_path / "manifest.csv"
    out_path = tmp_path / "frozen.json"
    table.to_csv(table_path, index=False)
    manifest.to_csv(manifest_path, index=False)

    script = SCRIPTS / "fit_frozen_gwtc_unbinned_kde_mode.py"
    subprocess.run(
        [
            sys.executable,
            str(script),
            "--table",
            str(table_path),
            "--manifest",
            str(manifest_path),
            "--output",
            str(out_path),
            "--bandwidth-multipliers",
            "0.8,1.0",
            "--k-min",
            "0.5",
            "--k-max",
            "5.0",
            "--n-k",
            "8",
            "--gh-n",
            "20",
        ],
        check=True,
    )

    frozen = json.loads(out_path.read_text(encoding="utf-8"))
    assert frozen["training_n_positive_finite"] == train_n
    assert frozen["train_ell_max"] < 3.0
    assert frozen["holdout_evaluated"] is False
