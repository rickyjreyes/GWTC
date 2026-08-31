import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from generate_gwtc4_population_null_sensitivity import (  # noqa: E402
    make_label_nulls,
    make_mass_precision_stratified_nulls,
    make_mass_residual_nulls,
    make_mass_stratified_nulls,
    quantile_codes,
    resolve_dynamic_subset,
)


def _events(n=24):
    mass = np.linspace(2.0, 5.0, n)
    spin = 0.45 + 0.08 * (mass - mass.mean()) + 0.03 * np.sin(4.0 * mass)
    width = np.linspace(0.04, 0.14, n)
    return pd.DataFrame(
        {
            "event": [f"GW{i:06d}_000000" for i in range(n)],
            "log_total_mass_source_median": mass,
            "final_spin_median": spin,
            "final_spin_width_16_84": width,
        }
    )


def test_label_null_preserves_global_spin_multiset():
    df = _events()
    y_perm, idx, meta = make_label_nulls(df, outer_n=5, seed=17)
    original = np.sort(df.loc[idx, "final_spin_median"].to_numpy(float))
    for row in y_perm:
        assert np.allclose(np.sort(row[idx]), original)
    assert meta["null_type"] == "final_spin_event_label_permutation"


def test_mass_stratified_null_preserves_each_mass_stratum_multiset():
    df = _events()
    y_perm, idx, _ = make_mass_stratified_nulls(
        df, outer_n=4, seed=19, mass_var="log_total_mass_source", mass_bins=4
    )
    mass = df["log_total_mass_source_median"].to_numpy(float)
    codes, _ = quantile_codes(mass, idx, 4)
    y = df["final_spin_median"].to_numpy(float)
    for row in y_perm:
        for code in np.unique(codes):
            local = idx[codes == code]
            assert np.allclose(np.sort(row[local]), np.sort(y[local]))


def test_mass_precision_stratified_preserves_global_multiset_and_reports_cells():
    df = _events()
    y_perm, idx, meta = make_mass_precision_stratified_nulls(
        df,
        outer_n=3,
        seed=23,
        mass_var="log_total_mass_source",
        mass_bins=3,
        precision_bins=2,
    )
    original = np.sort(df.loc[idx, "final_spin_median"].to_numpy(float))
    for row in y_perm:
        assert np.allclose(np.sort(row[idx]), original)
    assert meta["joint_cell_sizes"]
    assert sum(meta["joint_cell_sizes"].values()) == len(idx)


def test_mass_residual_null_is_deterministic_and_physically_bounded():
    df = _events()
    a, idx_a, meta_a = make_mass_residual_nulls(
        df,
        outer_n=5,
        seed=29,
        mass_var="log_total_mass_source",
        residual_degree=2,
    )
    b, idx_b, _ = make_mass_residual_nulls(
        df,
        outer_n=5,
        seed=29,
        mass_var="log_total_mass_source",
        residual_degree=2,
    )
    assert np.array_equal(idx_a, idx_b)
    assert np.allclose(a, b)
    assert np.all(np.abs(a[:, idx_a]) < 1.0)
    assert meta_a["generated_final_spin_min"] > -1.0
    assert meta_a["generated_final_spin_max"] < 1.0


def test_dynamic_quantile_subset_resolves_when_threshold_text_changes():
    expected = pd.DataFrame({"event": ["GW1"]})
    subset_map = {
        "high_final_spin_final_spin_ge_q0.75_0.712345": expected,
        "rank_top32_final_spin": expected,
    }
    got = resolve_dynamic_subset(
        subset_map, "high_final_spin_final_spin_ge_q0.75_0.800000"
    )
    assert got is expected
    assert resolve_dynamic_subset(subset_map, "rank_top32_final_spin") is expected
