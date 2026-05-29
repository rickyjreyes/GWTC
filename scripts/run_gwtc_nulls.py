#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_gwtc_nulls.py
=================

Generate and store the parametric Poisson bootstrap null distribution for one
declared primary cell (variable + subset + bin count). This is the primary
null used to compute p_global. PASS-grade claims require null_n >= 1000.

Outputs:
    tables/gwtc_nulls.csv   (one row per null replicate: T_null, plus the
                             observed T_obs and the computed global_p repeated
                             for provenance)
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
    ap = argparse.ArgumentParser(description="Parametric Poisson bootstrap null for one primary GWTC cell.")
    ap.add_argument("--table", default=os.path.join("tables", "gwtc_events_clean.csv"))
    ap.add_argument("--variable", default="M_chirp")
    ap.add_argument("--catalog", default=None)
    ap.add_argument("--cumulative", action="store_true")
    ap.add_argument("--source-class", default=None)
    ap.add_argument("--p-astro-min", type=float, default=0.5)
    ap.add_argument("--far-max", type=float, default=None)
    ap.add_argument("--observing-run", default=None)
    ap.add_argument("--bins", type=int, default=20)
    ap.add_argument("--k-min", type=float, default=0.5)
    ap.add_argument("--k-max", type=float, default=40.0)
    ap.add_argument("--n-k", type=int, default=120)
    ap.add_argument("--null-n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--baseline-degree", type=int, default=3)
    ap.add_argument("--out", default=os.path.join("tables", "gwtc_nulls.csv"))
    args = ap.parse_args(argv)

    df = L.load_clean_table(args.table)
    sub = L.select_subset(
        df, catalog=args.catalog, cumulative=args.cumulative,
        source_class=args.source_class, p_astro_min=args.p_astro_min,
        far_max=args.far_max, observing_run=args.observing_run,
    )
    values = L.get_variable_values(sub, args.variable)
    if values.size < 20:
        print(f"[nulls] too few events ({values.size}); INCOMPLETE.")
        pd.DataFrame(columns=["replicate", "T_null", "T_obs", "global_p"]).to_csv(args.out, index=False)
        return 0

    field = L.build_density_field(values, bin_count=args.bins)
    k_grid = L.make_k_grid(args.k_min, args.k_max, args.n_k)
    scan = L.scan_k(field, k_grid, baseline_degree=args.baseline_degree)
    null = L.run_poisson_null(
        field, k_grid, scan.delta_d_star, null_n=args.null_n,
        seed=args.seed, baseline_degree=args.baseline_degree, refit_baseline=True,
    )

    out = pd.DataFrame({
        "replicate": np.arange(args.null_n),
        "T_null": null.t_null,
        "T_obs": scan.delta_d_star,
        "global_p": null.global_p,
    })
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[nulls] T_obs={scan.delta_d_star:.3f}  k_star={scan.k_star:.3f}  "
          f"n_star={scan.n_star:.3f}  global_p={null.global_p:.4f}  null_n={args.null_n}")
    print(f"[nulls] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
