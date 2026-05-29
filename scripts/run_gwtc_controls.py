#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_gwtc_controls.py
====================

Negative-control discipline. The primary candidate must NOT be beaten by
controls that contain no genuine log-domain structure. We run:

  1. smooth_resample : replace event values with samples drawn from the smooth
                       fitted baseline population (mu0 over ell). Structure
                       here would indicate the scan manufactures modes.
  2. uniform_ell     : draw events uniformly in ell over the active support.
  3. jitter          : perturb each event within its (asymmetric) posterior
                       uncertainty and re-scan. A binning artifact is sensitive
                       to jitter; a robust mode is not.
  4. bounded_chi_eff : run on |chi_eff|, a bounded NON-scale coordinate, and
                       verify it does not create a false primary PASS.

Each control is scanned with the same k-grid and a Poisson null, exactly like
the primary. A primary candidate is Class I only if controls do NOT beat it
(controls should have larger / comparable global_p, i.e. weaker significance).

Outputs:
    tables/gwtc_control_results.csv
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wct_gwtc_lib as L  # noqa: E402


def _sample_from_baseline(values, n_draw, bins, rng):
    """Draw n_draw values from the smooth fitted baseline population in ell."""
    field = L.build_density_field(values, bin_count=max(bins, 30))
    mu0 = L.fit_baseline(field, degree=3)
    p = mu0 / mu0.sum()
    # sample bin indices then jitter uniformly within each bin
    idx = rng.choice(len(field.centers), size=n_draw, p=p)
    width = field.edges[1] - field.edges[0]
    ell = field.centers[idx] + rng.uniform(-width / 2, width / 2, size=n_draw)
    return np.exp(ell)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Negative-control battery for the GWTC log-domain diagnostic.")
    ap.add_argument("--table", default=os.path.join("tables", "gwtc_events_clean.csv"))
    ap.add_argument("--variable", default="M_chirp")
    ap.add_argument("--catalog", default=None)
    ap.add_argument("--cumulative", action="store_true")
    ap.add_argument("--source-class", default=None)
    ap.add_argument("--p-astro-min", type=float, default=0.5)
    ap.add_argument("--bins", type=int, default=20)
    ap.add_argument("--k-min", type=float, default=0.5)
    ap.add_argument("--k-max", type=float, default=40.0)
    ap.add_argument("--n-k", type=int, default=120)
    ap.add_argument("--null-n", type=int, default=400)
    ap.add_argument("--n-control-realizations", type=int, default=5,
                    help="Number of random realizations for resample/uniform/jitter controls.")
    ap.add_argument("--seed", type=int, default=2024)
    ap.add_argument("--out", default=os.path.join("tables", "gwtc_control_results.csv"))
    args = ap.parse_args(argv)

    df = L.load_clean_table(args.table)
    sub = L.select_subset(df, catalog=args.catalog, cumulative=args.cumulative,
                          source_class=args.source_class, p_astro_min=args.p_astro_min)
    catalog_version = args.catalog or ("cumulative" if args.cumulative else "all")
    subset = f"{catalog_version}_pa{args.p_astro_min}"
    k_grid = L.make_k_grid(args.k_min, args.k_max, args.n_k)
    rng = np.random.default_rng(args.seed)

    values = L.get_variable_values(sub, args.variable)
    n_events = values.size
    rows = []

    # ---- smooth_resample (multiple realizations) ----
    for r in range(args.n_control_realizations):
        vv = _sample_from_baseline(values, n_events, args.bins, rng)
        row = L.analyze_cell(vv, variable=args.variable, catalog_version=catalog_version,
                             subset=subset, bin_count=args.bins, k_grid=k_grid,
                             null_n=args.null_n, seed=args.seed + r,
                             notes=f"control=smooth_resample;real={r}")
        rows.append(row)

    # ---- uniform_ell (multiple realizations) ----
    field0 = L.build_density_field(values, bin_count=args.bins)
    for r in range(args.n_control_realizations):
        ell = rng.uniform(field0.ell_min, field0.ell_max, size=n_events)
        row = L.analyze_cell(np.exp(ell), variable=args.variable, catalog_version=catalog_version,
                             subset=subset, bin_count=args.bins, k_grid=k_grid,
                             null_n=args.null_n, seed=args.seed + 100 + r,
                             notes=f"control=uniform_ell;real={r}")
        rows.append(row)

    # ---- jitter within posterior uncertainty (if uncertainty columns exist) ----
    unc_map = {"M_chirp": ("M_chirp_lower", "M_chirp_upper"),
               "M_total": ("M_total_lower", "M_total_upper"),
               "D_L": ("D_L_lower", "D_L_upper")}
    if args.variable in unc_map and all(c in sub.columns for c in unc_map[args.variable]):
        col = L.SCALE_VARIABLES[args.variable]
        lo_c, hi_c = unc_map[args.variable]
        base = sub[col].to_numpy(dtype=float)
        lo = np.abs(sub[lo_c].to_numpy(dtype=float))
        hi = np.abs(sub[hi_c].to_numpy(dtype=float))
        ok = np.isfinite(base) & (base > 0)
        for r in range(args.n_control_realizations):
            sigma = np.where(rng.standard_normal(base.shape) < 0, lo, hi)
            sigma = np.nan_to_num(sigma, nan=0.0)
            jit = base + rng.standard_normal(base.shape) * sigma
            vv = jit[ok & np.isfinite(jit) & (jit > 0)]
            row = L.analyze_cell(vv, variable=args.variable, catalog_version=catalog_version,
                                 subset=subset, bin_count=args.bins, k_grid=k_grid,
                                 null_n=args.null_n, seed=args.seed + 200 + r,
                                 notes=f"control=jitter;real={r}")
            rows.append(row)
    else:
        print(f"[ctrl] no uncertainty columns for {args.variable}; jitter control skipped.")

    # ---- bounded non-scale coordinate (|chi_eff|) ----
    vv = L.get_variable_values(sub, "abs_chi_eff")
    row = L.analyze_cell(vv, variable="abs_chi_eff", catalog_version=catalog_version,
                         subset=subset, bin_count=args.bins, k_grid=k_grid,
                         null_n=args.null_n, seed=args.seed + 300,
                         notes="control=bounded_chi_eff")
    rows.append(row)

    out = pd.DataFrame(rows, columns=L.CSV_COLUMNS)
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    out.to_csv(args.out, index=False)

    gp = pd.to_numeric(out["global_p"], errors="coerce")
    best = gp.min()
    print(f"[ctrl] ran {len(out)} controls; best (smallest) control global_p = {best}")
    print(f"[ctrl] wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
