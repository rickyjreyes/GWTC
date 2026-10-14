#!/usr/bin/env python3
"""Audit OFFICIAL posterior-field availability; do not estimate absent data (#13).

This is a metadata/coverage gate, not an exploratory scan. Official units and
frame conventions must be explicitly established before using each observable.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import h5py
import pandas as pd
from gwtc_posterior_uncertainty import file_sha

CANDIDATES={
    "E_rad":("radiated_energy","radiated_energy_source"),
    "peak_luminosity":("peak_luminosity",),
}
def audit_field_names(field_names):
    names=set(field_names)
    out={}
    for var,candidates in CANDIDATES.items():
        found=[f for f in candidates if f in names]
        out[var]={"status":"PRESENT_EXACT_NAME" if len(found)==1 else
                  ("AMBIGUOUS_MULTIPLE_NAMES" if found else "ABSENT"),
                  "exact_fields_found":found}
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--posterior-manifest",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    mf=pd.read_csv(args.posterior_manifest,keep_default_na=False)
    needed={"event_name","path","dataset_path","sha256","source_url","weight_path"}
    if set(mf.columns)!=needed or mf.event_name.duplicated().any():
        raise SystemExit("Use a complete #11 official posterior manifest with unique names")
    rows=[]
    for r in mf.itertuples(index=False):
        if not str(r.source_url).startswith("https://"):
            raise SystemExit("Official source URL required")
        if file_sha(r.path)!=r.sha256:raise SystemExit("SHA256 mismatch: "+r.path)
        with h5py.File(r.path,"r") as f:
            if r.dataset_path not in f:raise SystemExit("Posterior dataset missing: "+r.event_name)
            names=f[r.dataset_path].dtype.names
        if names is None:raise SystemExit("Expected official structured posterior table")
        fields=audit_field_names(names)
        for observable,info in fields.items():
            rows.append({"event_name":r.event_name,"observable":observable,
                "status":info["status"],"field_names":";".join(info["exact_fields_found"]),
                "declared_dataset":r.dataset_path,"source_url":r.source_url,"file_sha256":r.sha256,
                "units_verified":False,"frame_convention_verified":False,
                "analysis_ready":False})
    if args.output.exists():raise SystemExit("Refusing overwrite")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output,index=False)
    print("Field-coverage audit only; units/frame/model provenance MUST be independently verified.")
if __name__=="__main__":main()
