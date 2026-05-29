#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_gwtc_threshold_stress.py
============================

Catalog-threshold and observing-run stress test. Repeat the scan for the same
variable under different selection branches:
    - p_astro >= 0.5
    - p_astro >= 0.9
    - FAR < 1/year (if FAR available)
    - BBH-only subset
    - per observing run (O1..O4b) and cumulative

A real universal structure should not appear only because one threshold or one
catalog update was chosen after the fact. Each branch is recorded so the search
family is auditable (look-elsewhere discipline).

Outputs:
    tables/gwtc_threshold_stress.csv
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wct_gwtc_lib as L  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Threshold / observing-run stress test.")
    ap.add_argument("--table", default=os.path.join("tables", "gwtc_events_clean.csv"))
    ap.add_argument("--variable", default="M_chirp")
    ap.add_argument("--bins", type=int, default=20)
    ap.add_argument("--k-min", type=float, default=0.5)
    ap.add_argument("--k-max", type=float, default=40.0)
    ap.add_argument("--n-k", type=int, default=120)
    ap.add_argument("--null-n", type=int, default=400)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", default=os.path.join("tables", "gwtc_threshold_stress.csv"))
    args = ap.parse_args(argv)

    df = L.load_clean_table(args.table)
    k_grid = L.make_k_grid(args.k_min, args.k_max, args.n_k)

    # Declared search branches (threshold + observing-run family).
    branches = [
        dict(label="cumulative_pa0.5", cumulative=True, p_astro_min=0.5),
        dict(label="cumulative_pa0.9", cumulative=True, p_astro_min=0.9),
        dict(label="cumulative_far1yr", cumulative=True, p_astro_min=0.0, far_max=1.0),
        dict(label="cumulative_BBH_pa0.5", cumulative=True, p_astro_min=0.5, source_class="BBH"),
        dict(label="O3a_pa0.5", observing_run="O3a", p_astro_min=0.5),
        dict(label="O3b_pa0.5", observing_run="O3b", p_astro_min=0.5),
        dict(label="O4a_pa0.5", observing_run="O4a", p_astro_min=0.5),
        dict(label="O4b_pa0.5", observing_run="O4b", p_astro_min=0.5),
    ]

    rows = []
    for br in branches:
        label = br.pop("label")
        sub = L.select_subset(df, **br)
        values = L.get_variable_values(sub, args.variable)
        catalog_version = "cumulative" if br.get("cumulative") else (br.get("observing_run") or "all")
        row = L.analyze_cell(
            values, variable=args.variable, catalog_version=catalog_version,
            subset=label, bin_count=args.bins, k_grid=k_grid, null_n=args.null_n,
            seed=args.seed, notes="threshold_stress",
        )
        rows.append(row)
        print(f"[thr] {label:22s} n_star={row['n_star']}  global_p={row['global_p']}  "
              f"verdict={row['verdict_label']}")

    out = pd.DataFrame(rows, columns=L.CSV_COLUMNS)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[thr] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
