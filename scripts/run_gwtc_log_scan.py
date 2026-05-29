#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_gwtc_log_scan.py
====================

Run the primary log-domain residual scan on one admissible scale variable and
one subset of the GWTC clean table. Writes the per-k scan curve and a canonical
one-row summary.

This is a WCT-motivated diagnostic only. It does not prove WCT.

Outputs:
    tables/gwtc_scan_primary.csv          (canonical one-row summary, appended)
    outputs/summary/scan_<var>_<subset>_B<bins>.csv   (per-k DeltaD curve)
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wct_gwtc_lib as L  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Primary GWTC log-domain residual scan on one variable/subset.")
    ap.add_argument("--table", default=os.path.join("tables", "gwtc_events_clean.csv"))
    ap.add_argument("--variable", default="M_chirp",
                    help="Scale variable: " + ", ".join(L.SCALE_VARIABLES) +
                         " (bounded diagnostic: " + ", ".join(L.BOUNDED_VARIABLES) + ").")
    ap.add_argument("--catalog", default=None, help="catalog_version filter, e.g. GWTC-5.0.")
    ap.add_argument("--cumulative", action="store_true", help="Use deduplicated primary entries.")
    ap.add_argument("--source-class", default=None, help="BBH | NSBH | BNS.")
    ap.add_argument("--p-astro-min", type=float, default=0.5)
    ap.add_argument("--far-max", type=float, default=None)
    ap.add_argument("--observing-run", default=None)
    ap.add_argument("--bins", type=int, default=20)
    ap.add_argument("--k-min", type=float, default=0.5)
    ap.add_argument("--k-max", type=float, default=40.0)
    ap.add_argument("--n-k", type=int, default=200)
    ap.add_argument("--null-n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--baseline-degree", type=int, default=3)
    ap.add_argument("--subset-label", default=None, help="Override subset label in output.")
    ap.add_argument("--out-summary", default=os.path.join("tables", "gwtc_scan_primary.csv"))
    ap.add_argument("--curve-dir", default=os.path.join("outputs", "summary"))
    ap.add_argument("--append", action="store_true", help="Append to summary CSV instead of overwriting.")
    args = ap.parse_args(argv)

    df = L.load_clean_table(args.table)
    sub = L.select_subset(
        df, catalog=args.catalog, cumulative=args.cumulative,
        source_class=args.source_class, p_astro_min=args.p_astro_min,
        far_max=args.far_max, observing_run=args.observing_run,
    )
    values = L.get_variable_values(sub, args.variable)

    subset_label = args.subset_label or _subset_label(args)
    catalog_version = args.catalog or ("cumulative" if args.cumulative else "all")
    k_grid = L.make_k_grid(args.k_min, args.k_max, args.n_k)

    row = L.analyze_cell(
        values, variable=args.variable, catalog_version=catalog_version,
        subset=subset_label, bin_count=args.bins, k_grid=k_grid,
        null_n=args.null_n, seed=args.seed, baseline_degree=args.baseline_degree,
        notes="primary_scan",
    )

    # Also write the per-k curve when there were enough events.
    if np.isfinite(row["DeltaD_star"]):
        field = L.build_density_field(values, bin_count=args.bins)
        scan = L.scan_k(field, k_grid, baseline_degree=args.baseline_degree)
        os.makedirs(args.curve_dir, exist_ok=True)
        curve_path = os.path.join(
            args.curve_dir, f"scan_{args.variable}_{subset_label}_B{args.bins}.csv")
        pd.DataFrame({
            "k": scan.k_grid,
            "DeltaD": scan.delta_d,
            "n": scan.k_grid * field.delta_ell_active / (2 * np.pi),
        }).to_csv(curve_path, index=False)
        print(f"[scan] wrote curve {curve_path}")

    _write_row(args.out_summary, row, append=args.append)
    print("[scan] result:")
    for k in L.CSV_COLUMNS:
        print(f"    {k:16s} {row[k]}")
    return 0


def _subset_label(args) -> str:
    parts = []
    if args.cumulative:
        parts.append("cumulative")
    if args.catalog:
        parts.append(args.catalog)
    if args.source_class:
        parts.append(args.source_class)
    if args.observing_run:
        parts.append(args.observing_run)
    parts.append(f"pa{args.p_astro_min}")
    if args.far_max is not None:
        parts.append(f"far{args.far_max}")
    return "_".join(parts) if parts else "all"


def _write_row(path, row, append=False):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    df_new = pd.DataFrame([row], columns=L.CSV_COLUMNS)
    if append and os.path.exists(path):
        old = pd.read_csv(path)
        df_new = pd.concat([old, df_new], ignore_index=True)
    df_new.to_csv(path, index=False)


if __name__ == "__main__":
    raise SystemExit(main())
