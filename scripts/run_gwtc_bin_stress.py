#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_gwtc_bin_stress.py
======================

Bin-count stress test. Re-run the log-domain scan on the same variable/subset
across several bin counts (default 20, 30, 40, 50). A candidate is stronger if
n_star (the active-domain winding number), not necessarily raw k_star, is
stable across bin counts. Aliasing-driven peaks above the per-binning Nyquist
frequency will move and produce unstable n_star -> fragile (Class III).

Outputs:
    tables/gwtc_bin_stress.csv
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
    ap = argparse.ArgumentParser(description="Bin-count stress test for one GWTC variable/subset.")
    ap.add_argument("--table", default=os.path.join("tables", "gwtc_events_clean.csv"))
    ap.add_argument("--variable", default="M_chirp")
    ap.add_argument("--catalog", default=None)
    ap.add_argument("--cumulative", action="store_true")
    ap.add_argument("--source-class", default=None)
    ap.add_argument("--p-astro-min", type=float, default=0.5)
    ap.add_argument("--bins", type=int, nargs="+", default=[20, 30, 40, 50])
    ap.add_argument("--k-min", type=float, default=0.5)
    ap.add_argument("--k-max", type=float, default=40.0)
    ap.add_argument("--n-k", type=int, default=120)
    ap.add_argument("--null-n", type=int, default=400)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--out", default=os.path.join("tables", "gwtc_bin_stress.csv"))
    args = ap.parse_args(argv)

    df = L.load_clean_table(args.table)
    sub = L.select_subset(df, catalog=args.catalog, cumulative=args.cumulative,
                          source_class=args.source_class, p_astro_min=args.p_astro_min)
    values = L.get_variable_values(sub, args.variable)
    catalog_version = args.catalog or ("cumulative" if args.cumulative else "all")
    subset = f"{catalog_version}_pa{args.p_astro_min}"
    if args.source_class:
        subset += f"_{args.source_class}"
    k_grid = L.make_k_grid(args.k_min, args.k_max, args.n_k)

    rows = []
    for B in args.bins:
        row = L.analyze_cell(
            values, variable=args.variable, catalog_version=catalog_version,
            subset=subset, bin_count=B, k_grid=k_grid, null_n=args.null_n,
            seed=args.seed, notes="bin_stress",
        )
        rows.append(row)
        print(f"[bin] B={B:3d}  k_best={row['k_best']}  n_star={row['n_star']}  "
              f"global_p={row['global_p']}  verdict={row['verdict_label']}")

    out = pd.DataFrame(rows, columns=L.CSV_COLUMNS)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.to_csv(args.out, index=False)

    n_stars = out["n_star"].to_numpy(dtype=float)
    n_stars = n_stars[np.isfinite(n_stars)]
    if n_stars.size >= 2:
        cv = np.std(n_stars) / np.mean(n_stars) if np.mean(n_stars) != 0 else float("inf")
        print(f"[bin] n_star spread: mean={np.mean(n_stars):.3f} "
              f"std={np.std(n_stars):.3f} CV={cv:.3f} "
              f"({'STABLE' if cv < 0.15 else 'UNSTABLE'})")
    print(f"[bin] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
