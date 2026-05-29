#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_gwtc_master_table.py
=========================

Concatenate every stage result CSV into one master table with a 'stage'
column, for auditing the full search family (look-elsewhere discipline).

Inputs (any that exist):
    tables/gwtc_scan_primary.csv
    tables/gwtc_bin_stress.csv
    tables/gwtc_variable_stress.csv
    tables/gwtc_threshold_stress.csv
    tables/gwtc_control_results.csv

Output:
    tables/gwtc_master_results.csv
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wct_gwtc_lib as L  # noqa: E402

STAGE_FILES = [
    ("primary", os.path.join("tables", "gwtc_scan_primary.csv")),
    ("bin_stress", os.path.join("tables", "gwtc_bin_stress.csv")),
    ("variable_stress", os.path.join("tables", "gwtc_variable_stress.csv")),
    ("threshold_stress", os.path.join("tables", "gwtc_threshold_stress.csv")),
    ("controls", os.path.join("tables", "gwtc_control_results.csv")),
]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Assemble the GWTC master results table.")
    ap.add_argument("--out", default=os.path.join("tables", "gwtc_master_results.csv"))
    args = ap.parse_args(argv)

    frames = []
    for stage, path in STAGE_FILES:
        if os.path.exists(path) and os.path.getsize(path) > 0:
            d = pd.read_csv(path)
            if d.empty:
                continue
            d.insert(0, "stage", stage)
            frames.append(d)
            print(f"[master] {stage}: {len(d)} rows from {path}")
        else:
            print(f"[master] {stage}: MISSING ({path})")

    if not frames:
        print("[master] no stage files found; nothing to assemble.")
        cols = ["stage"] + L.CSV_COLUMNS
        pd.DataFrame(columns=cols).to_csv(args.out, index=False)
        return 0

    master = pd.concat(frames, ignore_index=True)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    master.to_csv(args.out, index=False)
    print(f"[master] wrote {args.out}: {len(master)} total rows across {len(frames)} stages")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
