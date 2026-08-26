from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_baseline_v2_cv.py"


def test_baseline_v2_cli_writes_selected_training_model(tmp_path: Path) -> None:
    rng = np.random.default_rng(123)
    src = tmp_path / "events.csv"
    manifest = tmp_path / "manifest.csv"
    out = tmp_path / "cv.csv"

    ell = rng.normal(loc=2.5, scale=0.35, size=600)
    names = [f"E{i:04d}" for i in range(len(ell))]
    pd.DataFrame(
        {
            "commonName": names,
            "M_chirp": np.exp(ell),
            "p_astro": np.full(len(ell), 0.99),
            "is_primary_entry": np.ones(len(ell), dtype=bool),
        }
    ).to_csv(src, index=False)
    pd.DataFrame(
        {
            "event_name": names,
            "split": ["train"] * len(names),
        }
    ).to_csv(manifest, index=False)

    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--input",
            str(src),
            "--manifest",
            str(manifest),
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
    assert set(means["data_split"]) == {"train"}
    assert set(means["n_events"].astype(int)) == {600}
