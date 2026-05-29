#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_gwtc_variable_stress.py
===========================

Variable stress test. Run the same scan on each available admissible scale
variable: ln(M_chirp), ln(M_total), ln(D_L), ln(M_final), ln(redshift).
Bounded spin variables (|chi_eff|) are run separately and labelled as
diagnostic only; they can never produce a primary PASS.

Outputs:
    tables/gwtc_variable_stress.csv
"""

from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wct_gwtc_lib as L  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Variable stress test across admissible scale variables.")
    ap.add_argument("--table", default=os.path.join("tables", "gwtc_events_clean.csv"))
    ap.add_argument("--catalog", default=None)
    ap.add_argument("--cumulative", action="store_true")
    ap.add_argument("--source-class", default=None)
    ap.add_argument("--p-astro-min", type=float, default=0.5)
    ap.add_argument("--bins", type=int, default=20)
    ap.add_argument("--k-min", type=float, default=0.5)
    ap.add_argument("--k-max", type=float, default=40.0)
    ap.add_argument("--n-k", type=int, default=120)
    ap.add_argument("--null-n", type=int, default=400)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--include-bounded", action="store_true",
                    help="Also run bounded diagnostic variables (|chi_eff|).")
    ap.add_argument("--out", default=os.path.join("tables", "gwtc_variable_stress.csv"))
    args = ap.parse_args(argv)

    df = L.load_clean_table(args.table)
    sub = L.select_subset(df, catalog=args.catalog, cumulative=args.cumulative,
                          source_class=args.source_class, p_astro_min=args.p_astro_min)
    catalog_version = args.catalog or ("cumulative" if args.cumulative else "all")
    subset = f"{catalog_version}_pa{args.p_astro_min}"
    if args.source_class:
        subset += f"_{args.source_class}"
    k_grid = L.make_k_grid(args.k_min, args.k_max, args.n_k)

    variables = list(L.SCALE_VARIABLES)
    if args.include_bounded:
        variables += list(L.BOUNDED_VARIABLES)

    rows = []
    for var in variables:
        values = L.get_variable_values(sub, var)
        row = L.analyze_cell(
            values, variable=var, catalog_version=catalog_version, subset=subset,
            bin_count=args.bins, k_grid=k_grid, null_n=args.null_n,
            seed=args.seed, notes="variable_stress",
        )
        rows.append(row)
        print(f"[var] {var:12s} n_star={row['n_star']}  global_p={row['global_p']}  "
              f"verdict={row['verdict_label']}")

    out = pd.DataFrame(rows, columns=L.CSV_COLUMNS)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[var] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
