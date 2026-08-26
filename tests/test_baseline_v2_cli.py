from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_baseline_v2_cv.py"


def test_baseline_v2_cli_writes_selected_model(tmp_path: Path) -> None:
    rng = np.random.default_rng(123)
    src = tmp_path / "events.csv"
    out = tmp_path / "cv.csv"

    # Positive, smooth synthetic scale variable with enough events for 30 bins.
    ell = rng.normal(loc=2.5, scale=0.35, size=600)
    pd.DataFrame(
        {
            "M_chirp": np.exp(ell),
            "p_astro": np.full(len(ell), 0.99),
        }
    ).to_csv(src, index=False)

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--input",
            str(src),
            "--output",
            str(out),
            "--variable",
            "M_chirp",
            "--bins",
            "30",
            "--degrees",
            "3,4",
            "--folds",
            "5",
        ],
        check=True,
    )

    result = pd.read_csv(out)
    means = result[result["fold"].astype(str) == "mean"]
    assert len(means) == 2
    assert int(means["selected"].sum()) == 1
    assert set(means["degree"].astype(int)) == {3, 4}
