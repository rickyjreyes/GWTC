from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "make_holdout_manifest.py"


def test_holdout_manifest_is_frozen_by_catalog_prefix(tmp_path: Path) -> None:
    src = tmp_path / "events.csv"
    out = tmp_path / "manifest.csv"
    pd.DataFrame(
        {
            "commonName": ["A", "B", "C", "D"],
            "catalog_version": ["GWTC-4.1", "GWTC-5.0", "GWTC-4.0", "GWTC-5.0"],
            "p_astro": [0.99, 0.95, 0.7, 0.2],
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
            "--holdout-prefix",
            "GWTC-5",
            "--p-astro-min",
            "0.5",
        ],
        check=True,
    )

    got = pd.read_csv(out).set_index("event_name")
    assert got.loc["A", "split"] == "train"
    assert got.loc["B", "split"] == "holdout"
    assert got.loc["C", "split"] == "train"
    assert "D" not in got.index
    assert set(got["holdout_rule"]) == {"catalog_prefix:GWTC-5"}


def test_holdout_manifest_refuses_empty_holdout(tmp_path: Path) -> None:
    src = tmp_path / "events.csv"
    out = tmp_path / "manifest.csv"
    pd.DataFrame(
        {
            "commonName": ["A", "B"],
            "catalog_version": ["GWTC-4.0", "GWTC-4.1"],
            "p_astro": [0.9, 0.8],
        }
    ).to_csv(src, index=False)

    proc = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--input",
            str(src),
            "--output",
            str(out),
            "--holdout-prefix",
            "GWTC-5",
        ],
        text=True,
        capture_output=True,
    )
    assert proc.returncode != 0
    assert not out.exists()
