#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
make_gwtc_verdict.py
====================

Adjudicate the final harness verdict (PASS / PARTIAL / FAIL / INCOMPLETE) and
reliability class (I / II / III / IV) for the DECLARED primary cell, using the
full stress-test matrix. This is the only place PASS may be declared, and it is
declared from the global Poisson-bootstrap p-value plus robustness, never from
a local chi^2 p-value.

Verdict rules (declared a priori):

PASS  (Class I) requires ALL of:
    - declared primary variable is physically scale-like
    - primary global_p < 0.05
    - primary null_n >= 1000
    - bin stress shows stable n_star (coefficient of variation < CV_MAX)
    - threshold stress does not destroy the signal (>= 1 other threshold/run
      branch also significant, OR primary remains significant under p_astro>=0.9)
    - controls do NOT beat the primary (min control global_p >= primary global_p)
    - no arbitrary coordinate transform required (primary uses ln of a scale var)

FAIL  (Class IV-leaning) if:
    - primary global_p >= 0.05, OR
    - controls beat the primary, OR
    - the signal appears only under a single arbitrary branch

PARTIAL (Class II/III) if:
    - primary is significant but fails at least one robustness criterion

INCOMPLETE if:
    - data / nulls / provenance missing, or required stages not run

Output:
    tables/gwtc_verdict.csv
    outputs/summary/VERDICT.txt
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wct_gwtc_lib as L  # noqa: E402

CV_MAX = 0.15  # n_star coefficient-of-variation threshold for "stable".


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Adjudicate the final GWTC harness verdict and reliability class.")
    ap.add_argument("--primary-variable", default="M_chirp")
    ap.add_argument("--primary-subset", default="cumulative_pa0.5")
    ap.add_argument("--scan", default=os.path.join("tables", "gwtc_scan_primary.csv"))
    ap.add_argument("--bin-stress", default=os.path.join("tables", "gwtc_bin_stress.csv"))
    ap.add_argument("--threshold-stress", default=os.path.join("tables", "gwtc_threshold_stress.csv"))
    ap.add_argument("--controls", default=os.path.join("tables", "gwtc_control_results.csv"))
    ap.add_argument("--manifest", default=os.path.join("data", "manifest.csv"))
    ap.add_argument("--out", default=os.path.join("tables", "gwtc_verdict.csv"))
    ap.add_argument("--txt", default=os.path.join("outputs", "summary", "VERDICT.txt"))
    args = ap.parse_args(argv)

    reasons = []
    incomplete = []

    # ---- provenance present? ----
    provenance_ok = os.path.exists(args.manifest) and os.path.getsize(args.manifest) > 0
    if not provenance_ok:
        incomplete.append("provenance manifest missing")

    # ---- primary scan ----
    primary = None
    if os.path.exists(args.scan) and os.path.getsize(args.scan) > 0:
        sdf = pd.read_csv(args.scan)
        m = sdf[(sdf["variable"] == args.primary_variable)]
        if "subset" in sdf.columns:
            m2 = m[m["subset"] == args.primary_subset]
            m = m2 if not m2.empty else m
        if not m.empty:
            primary = m.iloc[-1].to_dict()
    else:
        incomplete.append("primary scan missing")

    if primary is None:
        incomplete.append("declared primary cell not found in scan results")

    # ---- gather metrics ----
    prim_gp = _num(primary["global_p"]) if primary else np.nan
    prim_nulln = _num(primary["null_n"]) if primary else np.nan
    prim_nstar = _num(primary["n_star"]) if primary else np.nan
    is_scale = args.primary_variable in L.SCALE_VARIABLES

    # bin-stress n_star stability
    n_star_cv = np.nan
    if os.path.exists(args.bin_stress) and os.path.getsize(args.bin_stress) > 0:
        bdf = pd.read_csv(args.bin_stress)
        bdf = bdf[bdf["variable"] == args.primary_variable]
        ns = pd.to_numeric(bdf["n_star"], errors="coerce").dropna().to_numpy()
        if ns.size >= 2 and np.mean(ns) != 0:
            n_star_cv = float(np.std(ns) / np.mean(ns))
    else:
        incomplete.append("bin stress missing")

    # threshold robustness: any OTHER branch significant?
    other_branch_sig = False
    p09_sig = False
    if os.path.exists(args.threshold_stress) and os.path.getsize(args.threshold_stress) > 0:
        tdf = pd.read_csv(args.threshold_stress)
        tdf = tdf[tdf["variable"] == args.primary_variable]
        tdf["gp"] = pd.to_numeric(tdf["global_p"], errors="coerce")
        others = tdf[tdf["subset"] != args.primary_subset]
        other_branch_sig = bool((others["gp"] < 0.05).any())
        p09 = tdf[tdf["subset"].astype(str).str.contains("pa0.9")]
        p09_sig = bool((pd.to_numeric(p09["global_p"], errors="coerce") < 0.05).any())
    else:
        incomplete.append("threshold stress missing")

    # Controls must be categorized, because "small global_p" means different
    # things for different controls:
    #   - structureless (smooth_resample, uniform_ell): a small p here means
    #     the scan MANUFACTURES significance on data with no real structure.
    #     If any structureless control beats the primary -> FAIL.
    #   - jitter: re-scan after perturbing within posterior uncertainty. A
    #     small p means the candidate SURVIVES jitter (robustness positive),
    #     so it is NOT used in the "controls beat primary" test.
    #   - bounded_chi_eff: a bounded NON-scale coordinate. If it is itself
    #     significant, the detected mode is not specific to scale-like
    #     variables -> the primary cannot be Class I.
    ctrl_struct_best = np.nan   # smallest structureless-control global_p
    bounded_coord_sig = False
    if os.path.exists(args.controls) and os.path.getsize(args.controls) > 0:
        cdf = pd.read_csv(args.controls)
        cdf["gp"] = pd.to_numeric(cdf["global_p"], errors="coerce")
        notes = cdf["notes"].astype(str)
        struct = cdf[notes.str.contains("smooth_resample") | notes.str.contains("uniform_ell")]
        if not struct.empty:
            ctrl_struct_best = float(struct["gp"].min())
        bounded = cdf[notes.str.contains("bounded_chi_eff") | (cdf["variable"] == "abs_chi_eff")]
        if not bounded.empty:
            bounded_coord_sig = bool((bounded["gp"] < 0.05).any())
    else:
        incomplete.append("controls missing")

    # ---- decide ----
    verdict = "INCOMPLETE"
    rclass = "IV"
    if not incomplete and primary is not None and np.isfinite(prim_gp):
        # A structureless control "beats" the primary only if it reaches a
        # global_p as small as the primary (manufacturing artifact).
        controls_beat = np.isfinite(ctrl_struct_best) and (ctrl_struct_best <= prim_gp)
        nstar_stable = np.isfinite(n_star_cv) and (n_star_cv < CV_MAX)
        threshold_ok = other_branch_sig or p09_sig
        nulln_ok = prim_nulln >= 1000
        sig = prim_gp < 0.05
        coord_specific = not bounded_coord_sig  # mode specific to scale vars?

        if sig:
            reasons.append(f"primary global_p={prim_gp:.4g} < 0.05")
        else:
            reasons.append(f"primary global_p={prim_gp:.4g} >= 0.05")
        reasons.append(f"null_n={int(prim_nulln)} ({'ok' if nulln_ok else 'TOO FEW (<1000)'})")
        reasons.append(f"n_star CV={n_star_cv:.3f} ({'stable' if nstar_stable else 'UNSTABLE across bin counts'})")
        reasons.append(f"threshold robustness: {'ok' if threshold_ok else 'fails'} "
                       f"(other_branch_sig={other_branch_sig}, p_astro>=0.9 sig={p09_sig})")
        reasons.append(f"structureless controls best global_p={ctrl_struct_best:.4g} "
                       f"({'STRUCTURELESS CONTROLS BEAT PRIMARY' if controls_beat else 'do not beat primary'})")
        reasons.append(f"bounded non-scale coordinate significant: {bounded_coord_sig} "
                       f"({'coordinate NOT specific to scale vars' if bounded_coord_sig else 'coordinate-specific'})")
        reasons.append(f"scale-like variable: {is_scale}")

        if not sig or controls_beat:
            # Either not significant, or significance is manufactured.
            verdict, rclass = "FAIL", "IV"
        elif (sig and is_scale and nulln_ok and nstar_stable and threshold_ok
              and not controls_beat and coord_specific):
            verdict, rclass = "PASS", "I"
        else:
            # Significant but missing a robustness layer.
            verdict = "PARTIAL"
            # Class III if bin/threshold/coordinate fragile, else Class II.
            fragile = (not nstar_stable) or (not threshold_ok) or (not coord_specific)
            rclass = "III" if fragile else "II"
    else:
        reasons.append("INCOMPLETE: " + "; ".join(incomplete) if incomplete else "INCOMPLETE")
        verdict, rclass = "INCOMPLETE", "IV"

    row = {
        "catalog_version": primary["catalog_version"] if primary else "",
        "variable": args.primary_variable,
        "subset": args.primary_subset,
        "primary_global_p": prim_gp,
        "primary_null_n": prim_nulln,
        "primary_n_star": prim_nstar,
        "n_star_cv_bin_stress": n_star_cv,
        "structureless_control_best_global_p": ctrl_struct_best,
        "bounded_coord_significant": bounded_coord_sig,
        "threshold_other_branch_sig": other_branch_sig,
        "is_scale_like": is_scale,
        "verdict_label": verdict,
        "reliability_class": rclass,
        "reasons": " | ".join(reasons),
    }
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    pd.DataFrame([row]).to_csv(args.out, index=False)

    os.makedirs(os.path.dirname(args.txt) or ".", exist_ok=True)
    with open(args.txt, "w") as fh:
        fh.write("WCT-MOTIVATED GWTC LOG-DOMAIN RESIDUAL DIAGNOSTIC -- VERDICT\n")
        fh.write("=" * 62 + "\n\n")
        fh.write("This is a diagnostic harness. It does not prove WCT and does not\n")
        fh.write("replace LVK population inference.\n\n")
        fh.write(f"Declared primary variable : {args.primary_variable} (ln-transform)\n")
        fh.write(f"Declared primary subset   : {args.primary_subset}\n\n")
        fh.write(f"VERDICT          : {verdict}\n")
        fh.write(f"RELIABILITY CLASS: {rclass}\n\n")
        fh.write("Reasoning:\n")
        for r in reasons:
            fh.write(f"  - {r}\n")
    print(f"[verdict] {verdict} (Class {rclass})")
    for r in reasons:
        print(f"    - {r}")
    print(f"[verdict] wrote {args.out} and {args.txt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
