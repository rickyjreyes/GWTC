#!/usr/bin/env python3
"""Provenance-checked, event-level posterior sensitivity for frozen V4 (#11).

This is a POSTERIOR-DRAW DIAGNOSTIC, NOT a hierarchical likelihood or a p-value:
one posterior draw is sampled per event in each realization, with event-wise
normalization, never by pooling thousands of samples as independent events.
Full PE-prior and injection corrections belong to issue #8.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

def file_sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fd:
        for chunk in iter(lambda: fd.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

def chirp_from_record(data):
    fields = data.dtype.names
    if fields is None: raise ValueError("Expected compound posterior_samples dataset")
    if "chirp_mass_source" in fields:
        z = np.asarray(data["chirp_mass_source"], dtype=float)
        source = "chirp_mass_source"
    elif {"mass_1_source", "mass_2_source"}.issubset(fields):
        m1 = np.asarray(data["mass_1_source"], float)
        m2 = np.asarray(data["mass_2_source"], float)
        z = (m1*m2)**(3.0/5.0) / (m1+m2)**(1.0/5.0)
        source = "derived_from_source_frame_component_masses"
    else:
        raise ValueError("Missing source-frame mass fields: never substitute detector-frame chirp mass")
    return z, source

def extract(file_path, dataset_path, expected_sha256, weight_path=""):
    if file_sha(file_path).lower() != str(expected_sha256).lower():
        raise ValueError("SHA256 mismatch for " + str(file_path))
    with h5py.File(file_path, "r") as f:
        if dataset_path not in f: raise ValueError("Declared posterior dataset absent")
        z, source = chirp_from_record(f[dataset_path][()])
        if weight_path:
            if weight_path not in f: raise ValueError("Declared posterior weight dataset absent")
            weights = np.asarray(f[weight_path][()], float)
        else:
            weights = np.ones(len(z), float)
    if z.ndim != 1 or weights.shape != z.shape:
        raise ValueError("Posterior arrays/weights require matching 1-D shape")
    valid = np.isfinite(z) & (z>0) & np.isfinite(weights) & (weights>=0)
    z, weights = z[valid], weights[valid]
    if not len(z) or float(weights.sum()) <= 0:
        raise ValueError("No finite positive weighted source-frame chirp mass")
    weights /= weights.sum()
    return np.log(z), weights, source

def draw_statistic(events, k, a, b, Z, n_draws, seed):
    """Independent one-draw-per-event replicate; vectorized event-wise weights."""
    if Z <= 0 or n_draws <= 0: raise ValueError("Invalid normalizer/draw count")
    rng = np.random.default_rng(seed)
    sums = np.zeros(n_draws, float)
    for ell, weights in events:
        idx = rng.choice(len(ell), size=n_draws, replace=True, p=weights)
        x = ell[idx]
        sums += a * np.cos(k*x) + b * np.sin(k*x) - np.log(Z)
    return 2.0*sums

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--posterior-manifest", required=True, type=Path,
                   help="CSV: event_name,path,dataset_path,sha256,source_url,weight_path")
    p.add_argument("--frozen", default="tables/gwtc_v4_frozen_structured_null.json", type=Path)
    p.add_argument("--table", default="tables/gwtc_events_clean.csv", type=Path)
    p.add_argument("--holdout-manifest", default="tables/gwtc_holdout_manifest.csv", type=Path)
    p.add_argument("--draws", type=int, default=10000)
    p.add_argument("--seed", type=int, default=20261014)
    p.add_argument("--output", default="tables/gwtc_posterior_uncertainty_diagnostic.json", type=Path)
    args = p.parse_args()
    freeze = json.loads(args.frozen.read_text())
    if freeze["schema"] != "GWTC_V4_STRUCTURED_NULL_FIXED_PREDICTION_V1" or float(freeze["k_pred"]) != 9.7:
        raise SystemExit("Expected exact frozen V4 model")
    if file_sha(args.table) != freeze["clean_table_sha256"] or file_sha(args.holdout_manifest) != freeze["manifest_sha256"]:
        raise SystemExit("Changed input table or holdout manifest")
    if args.draws < 2: raise SystemExit("At least two posterior-draw realizations required")
    manifest = pd.read_csv(args.posterior_manifest, keep_default_na=False)
    required = ["event_name","path","dataset_path","sha256","source_url","weight_path"]
    if list(manifest.columns) != required or manifest.event_name.duplicated().any():
        raise SystemExit("Manifest requires unique event_name and exact documented column order")
    manifest = manifest.sort_values("event_name")
    held = pd.read_csv(args.holdout_manifest)
    expected = set(held.loc[held.split=="holdout","event_name"].astype(str))
    input_names = set(manifest.event_name.astype(str))
    if not input_names.issubset(expected): raise SystemExit("Posterior manifest contains non-holdout event")
    if not input_names: raise SystemExit("Empty posterior manifest")
    table = pd.read_csv(args.table)
    primary = table.is_primary_entry.astype(str).str.lower().isin(("true","1","yes"))
    subset = table.loc[primary & table.commonName.astype(str).isin(expected)].copy()
    subset = subset.loc[np.isfinite(pd.to_numeric(subset.M_chirp, errors="coerce")) &
                        (pd.to_numeric(subset.M_chirp, errors="coerce")>0)]
    point = subset.set_index("commonName").M_chirp
    if point.index.duplicated().any(): raise SystemExit("Duplicate canonical primary input")
    expected_positive=set(point.index.astype(str))
    if not input_names.issubset(expected_positive):
        raise SystemExit("Posterior input contains events without finite positive point-estimate counterpart")
    mode=freeze["fixed_mode"]; k=float(mode["k"]); a=float(mode["a_cos"]); b=float(mode["b_sin"]); Z=float(mode["normalizer_Z"])
    events=[]; rows=[]
    for row in manifest.itertuples(index=False):
        if not str(row.source_url).startswith("https://"): raise SystemExit("Official source URL required for "+row.event_name)
        ell, w, source = extract(row.path,row.dataset_path,row.sha256,row.weight_path)
        events.append((ell,w))
        ess=float(1.0/np.sum(w*w))
        rows.append({"event_name":row.event_name,"path":row.path,"sha256":row.sha256,
          "dataset_path":row.dataset_path,"weight_path":row.weight_path,"source_url":row.source_url,
          "source_frame_field":source,"posterior_n":len(ell),"weight_ess":ess,
          "posterior_chirp_median":float(np.quantile(np.exp(ell),.5)),
          "point_chirp":float(point.loc[row.event_name])})
    scores=draw_statistic(events,k,a,b,Z,args.draws,args.seed)
    ell_point=np.log(point.loc[manifest.event_name.astype(str)].to_numpy(float))
    t_point=float(2*np.sum(a*np.cos(k*ell_point)+b*np.sin(k*ell_point)-np.log(Z)))
    out={"schema":"GWTC_V4_POSTERIOR_SENSITIVITY_V1",
         "scientific_status":"POSTERIOR_DRAW_DIAGNOSTIC_NOT_HIERARCHICAL_NOT_CALIBRATED",
         "frozen_sha256":file_sha(args.frozen),"manifest_sha256":file_sha(args.posterior_manifest),
         "table_sha256":file_sha(args.table),"holdout_manifest_sha256":file_sha(args.holdout_manifest),
         "holdout_manifest_rows":len(expected),"holdout_expected_positive_mass":len(expected_positive),
         "input_events":len(events),
         "missing_event_names":sorted(expected_positive-input_names),
         "coverage":"COMPLETE" if input_names==expected_positive else "INCOMPLETE",
         "point_statistic_on_matched_events":t_point,
         "draws":args.draws,"seed":args.seed,
         "posterior_draw_score_mean":float(scores.mean()),
         "posterior_draw_score_quantiles_05_50_95":np.quantile(scores,[.05,.5,.95]).tolist(),
         "posterior_draw_fraction_statistic_positive":float(np.mean(scores>0)),
         "calibrated_p":None,
         "events":rows}
    if args.output.exists():raise SystemExit("Refusing overwrite; preserve first result and rerun under a new filename")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print("Wrote",args.output,"coverage=",out["coverage"],"posterior test is diagnostic, NOT a p-value")

if __name__=="__main__":
    main()
