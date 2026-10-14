"""#16 full-family calibration: local != global; correlation must be retained."""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from calibrate_exploratory_family import calibrate_family


def test_rejects_missing_and_extra_cells():
    obs = pd.DataFrame({"cell_id": ["a", "b"], "statistic": [3., 3.]})
    with pytest.raises(ValueError, match="EXACTLY"):
        calibrate_family(obs, pd.DataFrame({"a": [1., 2.]}))
    with pytest.raises(ValueError, match="EXACTLY"):
        calibrate_family(obs, pd.DataFrame({"a": [1., 2.], "b": [2., 3.], "c": [4., 5.]}))


def test_joint_null_retains_perfect_dependence():
    values = np.arange(1., 21.)
    obs = pd.DataFrame({"cell_id": ["a", "b"], "statistic": [19.5, 19.5]})
    _, both, nulls = calibrate_family(obs, pd.DataFrame({"a": values, "b": values}))
    _, one, _ = calibrate_family(obs.iloc[:1], pd.DataFrame({"a": values}))
    assert both["family_global_p"] == one["family_global_p"]
    assert nulls.min_local_p.between(0, 1).all()


def test_full_family_includes_late_cell():
    obs = pd.DataFrame({"cell_id": ["a", "b"], "statistic": [0., 11.]})
    matrix = pd.DataFrame({"a": np.arange(10.), "b": np.arange(10.)})
    rows, summary, _ = calibrate_family(obs, matrix)
    assert summary["best_cell_ids"] == ["b"]
    assert rows.loc[0, "local_p"] == 1.0
    assert rows.loc[1, "local_p"] == pytest.approx(1 / 11)
    assert summary["family_global_p"] >= summary["best_local_p"]


def test_ties_are_conservative():
    obs = pd.DataFrame({"cell_id": ["a"], "statistic": [1.]})
    _, _, nulls = calibrate_family(obs, pd.DataFrame({"a": [1., 1., 1., 1.]}))
    assert (nulls.min_local_p == 1).all()
