import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from generate_gwtc4_population_null_matrix import (  # noqa: E402
    is_dynamic_subset,
    leave_one_out_rank_p,
    prepare_scan,
    score_prepared_one,
    upper_rank_p_observed,
)
from gwtc4_wct_subset_scan_compound import harmonic_scan  # noqa: E402


def _synthetic_subset(n=24, seed=7):
    import pandas as pd

    rng = np.random.default_rng(seed)
    x = np.sort(rng.uniform(1.2, 4.8, n))
    y = 0.15 * x + 0.3 * np.cos(6.2 * x) + rng.normal(0.0, 0.08, n)
    width = rng.uniform(0.05, 0.15, n)
    return pd.DataFrame(
        {
            "event": [f"GW{i:06d}_000000" for i in range(n)],
            "log_total_mass_source_median": x,
            "final_spin_median": y,
            "final_spin_width_16_84": width,
        }
    )


def test_fast_scan_matches_historical_lstsq_loop():
    df = _synthetic_subset()
    k_grid = np.linspace(0.5, 20.0, 160)
    prep = prepare_scan(
        df,
        "log_total_mass_source",
        "final_spin",
        min_events=12,
        k_grid=k_grid,
        degree=2,
    )
    y = df.set_index("event").loc[prep.events, "final_spin_median"].to_numpy(float)
    k_fast, d_fast = score_prepared_one(prep, y)

    x = prep.x
    w = prep.w
    scan_df, best, _, _ = harmonic_scan(x, y, w, k_grid, degree=2)

    assert np.isclose(k_fast, best["k_best"], atol=0.0, rtol=0.0)
    assert np.isclose(d_fast, best["delta_chi2"], atol=1e-10, rtol=1e-10)
    assert np.isclose(d_fast, scan_df["delta_chi2"].max(), atol=1e-10, rtol=1e-10)


def test_dynamic_subset_identification():
    assert is_dynamic_subset("high_final_spin_final_spin_ge_q0.75_0.8")
    assert is_dynamic_subset("low_final_spin_final_spin_le_q0.25_0.5")
    assert is_dynamic_subset("rank_top32_final_spin")
    assert not is_dynamic_subset("rank_top32_eta_gap_from_equal_mass")
    assert not is_dynamic_subset("high_mass_total_mass_source_ge_q0.75_70")


def test_outer_observed_rank_uses_plus_one():
    null_delta = np.array([[1.0, 4.0], [2.0, 3.0], [3.0, 2.0], [4.0, 1.0]])
    obs = np.array([3.5, 4.5])
    p = upper_rank_p_observed(null_delta, obs)
    assert np.allclose(p, [2.0 / 5.0, 1.0 / 5.0])


def test_leave_one_out_null_ranks_are_exchangeable_upper_tails():
    null_delta = np.array([[4.0], [3.0], [2.0], [1.0]])
    p = leave_one_out_rank_p(null_delta)
    assert np.allclose(p[:, 0], [0.25, 0.50, 0.75, 1.00])
