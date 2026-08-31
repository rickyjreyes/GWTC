import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from run_gwtc_population_global_null import (  # noqa: E402
    HISTORICAL_EDGES,
    aggregate_summary,
    empirical_upper_tail,
    histogram_chi2,
    validate_pvalues,
)


def test_historical_edges_match_reported_expectations_for_168_scans():
    widths = np.diff(HISTORICAL_EDGES)
    expected = 168 * widths
    assert len(widths) == 11
    assert np.isclose(expected[0], 8.4)
    assert np.isclose(expected[1], 8.4)
    assert np.allclose(expected[2:], 16.8)


def test_uniform_midpoints_have_zero_histogram_chi2():
    # Put counts in each bin exactly proportional to its width for N=220:
    # 11 in each 0.05 bin and 22 in each 0.10 bin.
    values = []
    for left, right in zip(HISTORICAL_EDGES[:-1], HISTORICAL_EDGES[1:]):
        count = 11 if np.isclose(right - left, 0.05) else 22
        values.extend([0.5 * (left + right)] * count)
    values = np.asarray(values, dtype=float)
    stat, counts, expected = histogram_chi2(values)
    assert values.size == 220
    assert np.allclose(counts, expected)
    assert np.isclose(stat, 0.0)


def test_aggregate_summary_counts_thresholds():
    p = np.concatenate(
        [
            np.full(16, 0.01),
            np.full(31, 0.07),
            np.full(121, 0.50),
        ]
    )
    summary = aggregate_summary(p)
    assert summary["n_scans"] == 168
    assert summary["count_p_lt_0p05"] == 16
    assert summary["count_p_lt_0p10"] == 47
    assert summary["poisson_tail_p_lt_0p10"] < 1e-10


def test_empirical_tail_uses_plus_one_correction():
    null_stats = np.array([1.0, 2.0, 3.0, 4.0])
    ge, p = empirical_upper_tail(null_stats, 3.0)
    assert ge == 2
    assert np.isclose(p, 3.0 / 5.0)


def test_validate_pvalues_rejects_out_of_range():
    try:
        validate_pvalues(np.array([0.1, 1.2]), "test")
    except SystemExit:
        pass
    else:
        raise AssertionError("Expected out-of-range p-value rejection")
