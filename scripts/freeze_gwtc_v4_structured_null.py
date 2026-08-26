#!/usr/bin/env python3
"""Freeze the GWTC V4 fixed-k structured-null robustness challenge.

V4 is not a new discovery scan.  The prospective prediction k=9.7 was already
published, so this script contains no k-grid arguments and never inspects the
holdout.  It performs two training-only operations:

1. fit the residual coefficients a,b at exactly k=9.7 using the already frozen
   V3 KDE baseline;
2. fit several non-periodic broken-power-law + Gaussian-peak population models
   under predeclared selection-weight scenarios.

The resulting JSON must be hashed/committed before the V4 evaluator is run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gwtc_structured_null import fit_structured_null  # noqa: E402
from gwtc_unbinned_kde import fit_mode_at_k  # noqa: E402

K_PRED = 9.7
K_BAND = (9.5, 10.0)
DEFAULT_SELECTION_GAMMAS = (0.0, 1.5, 2.5)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def event_column(df: pd.DataFrame) -> str:
    for name in ("commonName", "common_name", "event", "event_name", "name"):
        if name in df.columns:
            return name
    raise SystemExit("No event-name column found in clean table.")


def primary_mask(df: pd.DataFrame) -> pd.Series:
    if "is_primary_entry" not in df.columns:
        raise SystemExit("Clean table lacks is_primary_entry; refusing noncanonical V4 freeze.")
    s = df["is_primary_entry"]
    if pd.api.types.is_bool_dtype(s):
        return s.fillna(False)
    return s.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--table", default="tables/gwtc_events_clean.csv")
    p.add_argument("--manifest", default="tables/gwtc_holdout_manifest.csv")
    p.add_argument("--v3-frozen", default="tables/gwtc_v3_frozen_unbinned_kde_mode.json")
    p.add_argument("--selection-gammas", default="0,1.5,2.5")
    p.add_argument("--gh-n", type=int, default=40)
    p.add_argument("--grid-n", type=int, default=2048)
    p.add_argument("--output", default="tables/gwtc_v4_frozen_structured_null.json")
    args = p.parse_args()

    table_path = Path(args.table)
    manifest_path = Path(args.manifest)
    v3_path = Path(args.v3_frozen)

    v3 = json.loads(v3_path.read_text(encoding="utf-8"))
    if v3.get("schema") != "GWTC_V3_UNBINNED_KDE_FIXED_MODE_V1":
        raise SystemExit("Unrecognized V3 frozen-model schema.")
    if sha256_file(manifest_path) != v3.get("manifest_sha256"):
        raise SystemExit("Manifest differs from the one used to freeze V3.")
    if sha256_file(table_path) != v3.get("clean_table_sha256"):
        raise SystemExit("Clean table differs from the one used to freeze V3.")

    df = pd.read_csv(table_path)
    manifest = pd.read_csv(manifest_path)
    if not {"event_name", "split"}.issubset(manifest.columns):
        raise SystemExit("Manifest must contain event_name and split columns.")
    if not (manifest["split"] == "train").any() or not (manifest["split"] == "holdout").any():
        raise SystemExit("Manifest must contain both train and holdout rows.")

    variable = str(v3["variable"])
    if variable not in df.columns:
        raise SystemExit(f"Variable {variable!r} not found in clean table.")

    train_names = set(manifest.loc[manifest["split"] == "train", "event_name"].astype(str))
    ev_col = event_column(df)
    canonical = df.loc[primary_mask(df)].copy()
    canonical[ev_col] = canonical[ev_col].astype(str)
    train = canonical[canonical[ev_col].isin(train_names)].copy()

    z = pd.to_numeric(train[variable], errors="coerce").to_numpy(float)
    z = z[np.isfinite(z) & (z > 0.0)]
    if z.size < 30:
        raise SystemExit(f"Too few positive finite training {variable} values: {z.size}")
    ell = np.log(z)

    v3_centers = np.asarray(v3["train_ell"], dtype=float)
    if ell.size != v3_centers.size or not np.allclose(
        np.sort(ell), np.sort(v3_centers), rtol=0.0, atol=1e-12
    ):
        raise SystemExit("Training values do not match the frozen V3 KDE centers.")

    bandwidth = float(v3["selected_bandwidth"])
    fixed_mode = fit_mode_at_k(
        ell,
        v3_centers,
        bandwidth,
        K_PRED,
        gh_n=args.gh_n,
    )

    gammas = [float(x.strip()) for x in args.selection_gammas.split(",") if x.strip()]
    if not gammas:
        raise SystemExit("At least one selection gamma is required.")
    if len(set(gammas)) != len(gammas):
        raise SystemExit("Selection gamma list contains duplicates.")
    if any(g < 0.0 or not np.isfinite(g) for g in gammas):
        raise SystemExit("Selection gammas must be finite and non-negative.")

    scenarios = []
    for gamma in gammas:
        params = fit_structured_null(z, gamma, grid_n=args.grid_n)
        scenarios.append(
            {
                "name": f"broken_powerlaw_plus_peak_gamma_{gamma:g}",
                "selection_gamma": float(gamma),
                "population_family": "continuous_broken_powerlaw_plus_truncated_gaussian_peak",
                "selection_family": "relative_efficiency_(Mchirp/Mchirp_max)^gamma",
                "params": params.to_dict(),
            }
        )

    delta_ell = float(np.max(ell) - np.min(ell))
    n_equiv = float(K_PRED * delta_ell / (2.0 * math.pi))
    mass_ratio = float(math.exp(2.0 * math.pi / K_PRED))

    payload = {
        "schema": "GWTC_V4_STRUCTURED_NULL_FIXED_PREDICTION_V1",
        "scientific_status": "ROBUSTNESS_CHALLENGE_NOT_INDEPENDENT_REPLICATION",
        "reason": "GWTC-5 was already inspected in V2/V3. V4 tests whether training-fitted non-periodic structured populations plus selection can mimic the already published fixed k=9.7 residual.",
        "variable": variable,
        "coordinate": "ell=log(variable)",
        "prediction_status": "EXTERNALLY_PUBLISHED_AND_FROZEN_BEFORE_V4",
        "k_pred": K_PRED,
        "k_acceptance_band": [K_BAND[0], K_BAND[1]],
        "equivalent_mass_ratio_exp_2pi_over_k": mass_ratio,
        "n_equiv_over_training_span": n_equiv,
        "manifest_path": str(manifest_path).replace("\\", "/"),
        "manifest_sha256": sha256_file(manifest_path),
        "clean_table_path": str(table_path).replace("\\", "/"),
        "clean_table_sha256": sha256_file(table_path),
        "v3_frozen_path": str(v3_path).replace("\\", "/"),
        "v3_frozen_sha256": sha256_file(v3_path),
        "training_manifest_rows": int((manifest["split"] == "train").sum()),
        "training_n_positive_finite": int(z.size),
        "training_z_min": float(np.min(z)),
        "training_z_max": float(np.max(z)),
        "training_ell_min": float(np.min(ell)),
        "training_ell_max": float(np.max(ell)),
        "kde_bandwidth_frozen_from_v3": bandwidth,
        "gh_n": int(args.gh_n),
        "structured_null_grid_n": int(args.grid_n),
        "fixed_mode": {
            "k": float(fixed_mode.k),
            "a_cos": float(fixed_mode.a),
            "b_sin": float(fixed_mode.b),
            "amplitude": float(fixed_mode.amplitude),
            "phase_atan2_b_a": float(fixed_mode.phase_atan2_b_a),
            "normalizer_Z": float(fixed_mode.normalizer),
            "train_delta_2logl": float(fixed_mode.train_delta_2logl),
        },
        "structured_null_scenarios": scenarios,
        "holdout_evaluated_by_this_script": False,
        "interpretation_limit": "The structured null is phenomenological and selection weighting is approximate; survival does not replace full hierarchical LVK population inference with injection-based selection functions and event posteriors.",
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"Wrote {out}")
    print(f"Training N       : {z.size}")
    print(f"Published k      : {K_PRED:.8g} (NO SCAN)")
    print(f"Frozen a,b       : {fixed_mode.a:.8g}, {fixed_mode.b:.8g}")
    print(f"Frozen amplitude : {fixed_mode.amplitude:.8g}")
    print(f"Train Delta2logL : {fixed_mode.train_delta_2logl:.8g}")
    print(f"Mass ratio       : {mass_ratio:.8g}")
    for scenario in scenarios:
        q = scenario["params"]
        print(
            "Null scenario     : "
            f"gamma={scenario['selection_gamma']:g}, "
            f"alpha=({q['alpha_low']:.4g},{q['alpha_high']:.4g}), "
            f"break={q['break_z']:.4g}, peak_frac={q['peak_fraction']:.4g}, "
            f"mu={q['peak_mu']:.4g}, sigma={q['peak_sigma']:.4g}"
        )
    print("HOLDOUT HAS NOT BEEN EVALUATED BY THIS SCRIPT.")
    print("Hash/commit this V4 JSON before running the V4 evaluator.")


if __name__ == "__main__":
    main()
