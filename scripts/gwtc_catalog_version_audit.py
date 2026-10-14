#!/usr/bin/env python3
"""Version-specific frozen-mode sensitivity, NOT fresh replication (#12).

GWOSC version rows may be incremental, not three interchangeable cumulative
population snapshots. Absence from a subsequent version is NOT a removal.
Retain raw event diffs and report a narrow fixed-KDE reference calibration only.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import sys
import numpy as np
import pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from gwtc_unbinned_kde import fixed_model_null, fixed_log_likelihood_ratio

VERSIONS=("GWTC-4.0","GWTC-4.1","GWTC-5.0")
FIELDS=("M_chirp","M_total","mass_1_source","mass_2_source","M_final","redshift","p_astro","far")
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def version_index(df,version):
    sub=df.loc[df.catalog_version==version].copy()
    if sub.commonName.isna().any() or sub.commonName.duplicated().any():
        raise ValueError("Ambiguous/missing event identity within version "+version)
    return sub.set_index("commonName",drop=False)
def version_diff(a,b,aa,bb):
    """Missing is 'not listed in this version' rather than physical removal."""
    rows=[]
    for name in sorted(set(a.index)|set(b.index)):
        if name not in a.index:
            rows.append({"version_a":aa,"version_b":bb,"event_name":name,"status":"not_listed_in_a","changed_fields":""});continue
        if name not in b.index:
            rows.append({"version_a":aa,"version_b":bb,"event_name":name,"status":"not_listed_in_b","changed_fields":""});continue
        changed=[]
        for fld in FIELDS:
            x=pd.to_numeric(pd.Series([a.loc[name,fld]]),errors="coerce").iloc[0]
            y=pd.to_numeric(pd.Series([b.loc[name,fld]]),errors="coerce").iloc[0]
            if pd.isna(x) and pd.isna(y):continue
            if pd.isna(x) != pd.isna(y) or not np.isclose(x,y,rtol=0,atol=1e-12):
                changed.append(fld)
        rows.append({"version_a":aa,"version_b":bb,"event_name":name,
                     "status":"revised" if changed else "unchanged","changed_fields":";".join(changed)})
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table",default="tables/gwtc_events_clean.csv",type=Path)
    p.add_argument("--v3",default="tables/gwtc_v3_frozen_unbinned_kde_mode.json",type=Path)
    p.add_argument("--v4",default="tables/gwtc_v4_frozen_structured_null.json",type=Path)
    p.add_argument("--null-n",type=int,default=10000)
    p.add_argument("--seed",type=int,default=20261014)
    p.add_argument("--out-prefix",default="tables/gwtc_catalog_version_audit",type=Path)
    args=p.parse_args()
    v4=json.loads(args.v4.read_text());v3=json.loads(args.v3.read_text())
    if sha(args.table)!=v4["clean_table_sha256"] or sha(args.v3)!=v4["v3_frozen_sha256"]:
        raise SystemExit("Input drift from frozen V4")
    if args.null_n < 1000:raise SystemExit("Require >=1000 nulls for version-reference diagnostics")
    df=pd.read_csv(args.table)
    versions={ver:version_index(df,ver) for ver in VERSIONS}
    diffs=[]
    for aa,bb in combinations(VERSIONS,2):
        diffs.extend(version_diff(versions[aa],versions[bb],aa,bb))
    m=v4["fixed_mode"]; out=[]
    for i,(version,sub) in enumerate(versions.items()):
        selected=sub.loc[(pd.to_numeric(sub.p_astro,errors="coerce")>=.5)].copy()
        z=pd.to_numeric(selected.M_chirp,errors="coerce").to_numpy(float)
        z=z[np.isfinite(z)&(z>0)]
        row={"version":version,"raw_version_entries":len(sub),"positive_pastro_and_mass_n":len(z),
             "zero_or_missing_mass_after_cut":len(selected)-len(z),
             "statistic_type":"frozen V4 k=9.7 score against V3 training KDE reference",
             "evidence_status":"DESCRIPTIVE_CROSS_VERSION_SENSITIVITY_NOT_PROSPECTIVE",
             "p_value_not_population_corrected":True}
        if len(z):
            ell=np.log(z)
            delta=fixed_log_likelihood_ratio(ell,m["k"],m["a_cos"],m["b_sin"],m["normalizer_Z"])
            null=fixed_model_null(kde_centers=v3["train_ell"],bandwidth=v3["selected_bandwidth"],
                    n_events=len(ell),k=m["k"],a=m["a_cos"],b=m["b_sin"],normalizer=m["normalizer_Z"],
                    null_n=args.null_n,seed=args.seed+i*1000003)
            row.update({"statistic_2loglr":delta,"null_n":args.null_n,"null_exceed":int((null>=delta).sum()),
                "reference_baseline_p":float((1+(null>=delta).sum())/(1+args.null_n))})
        out.append(row)
    args.out_prefix.parent.mkdir(parents=True,exist_ok=True)
    if Path(str(args.out_prefix)+"_results.json").exists():raise SystemExit("Refusing overwrite of version-audit result")
    Path(str(args.out_prefix)+"_results.json").write_text(json.dumps({
      "schema":"GWTC_VERSION_SENSITIVITY_V1","table_sha256":sha(args.table),
      "v3_sha256":sha(args.v3),"v4_sha256":sha(args.v4),"results":out,
      "limitation":"These catalog-version entries may be incremental and overlapping; reference KDE calibration is not a cross-version prospective replication or full hierarchical population test."
    },indent=2)+"\n")
    pd.DataFrame(diffs).to_csv(str(args.out_prefix)+"_event_diff.csv",index=False)
    print("Wrote",str(args.out_prefix)+"_results.json","and event_diff.csv")
if __name__=="__main__":main()
