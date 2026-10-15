#!/usr/bin/env python3
"""Exact additive decomposition of ARCHIVED tagged-version fixed-V4 score changes (#12).

This does not reconstruct official cumulative releases, establish prospective
independence, or assign a population-corrected p-value. Compare physical event
identities without treating not-listed as deletion. Only frozen V4 coefficients,
frozen V3 provenance and the original clean-table SHA are accepted.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from gwtc_catalog_version_audit import VERSIONS,version_index
from gwtc_source_provenance import validate_frozen_csv

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def scored_versions(table,mode):
    versions={}
    for ver in VERSIONS:
        sub=version_index(table,ver)
        z=pd.to_numeric(sub.M_chirp,errors="coerce")
        p=pd.to_numeric(sub.p_astro,errors="coerce")
        selected=(p>=.5)&np.isfinite(z)&(z>0)
        contribution=pd.Series(0.,index=sub.index,dtype=float)
        ell=np.log(z.loc[selected].to_numpy(float))
        contribution.loc[selected]=2*(float(mode["a_cos"])*np.cos(float(mode["k"])*ell)
            +float(mode["b_sin"])*np.sin(float(mode["k"])*ell)
            -np.log(float(mode["normalizer_Z"])))
        versions[ver]={"source":sub,"selected":selected,"score":contribution}
    return versions

def decompose(versions):
    """Per-ID ledger, with exact identity matching; no fuzzy rename or non-GWTC proxy."""
    ledger=[];summaries=[]
    for va,vb in combinations(VERSIONS,2):
        a,b=versions[va],versions[vb]
        keys=sorted(set(a["source"].index)|set(b["source"].index))
        categories={"shared_score_change":0.,"not_listed_in_a":0.,"not_listed_in_b":0.}
        for key in keys:
            ina=key in a["source"].index;inb=key in b["source"].index
            sa=float(a["score"].loc[key]) if ina else 0.
            sb=float(b["score"].loc[key]) if inb else 0.
            selected_a=bool(a["selected"].loc[key]) if ina else False
            selected_b=bool(b["selected"].loc[key]) if inb else False
            category=("shared_score_change" if ina and inb else
                      "not_listed_in_a" if inb else "not_listed_in_b")
            categories[category]+=sb-sa
            ledger.append({"version_a":va,"version_b":vb,"commonName":key,
                "identity_category":category,"selected_a":selected_a,"selected_b":selected_b,
                "mass_a":float(a["source"].loc[key,"M_chirp"]) if ina and pd.notna(a["source"].loc[key,"M_chirp"]) else None,
                "mass_b":float(b["source"].loc[key,"M_chirp"]) if inb and pd.notna(b["source"].loc[key,"M_chirp"]) else None,
                "contribution_a":sa,"contribution_b":sb,"score_change_b_minus_a":sb-sa,
                "selection_flip":ina and inb and selected_a!=selected_b})
        score_a=float(a["score"].sum());score_b=float(b["score"].sum())
        delta=score_b-score_a
        if not np.isclose(sum(categories.values()),delta,atol=1e-10,rtol=0):
            raise ValueError("Category decomposition fails score conservation")
        summaries.append({"version_a":va,"version_b":vb,"selected_a":int(a["selected"].sum()),
            "selected_b":int(b["selected"].sum()),"score_a":score_a,"score_b":score_b,
            "score_change_b_minus_a":delta,"decomposition":categories,
            "common_id_count":len(set(a["source"].index)&set(b["source"].index)),
            "unmatched_id_a":len(set(a["source"].index)-set(b["source"].index)),
            "unmatched_id_b":len(set(b["source"].index)-set(a["source"].index)),
            "shared_selection_flips":sum(bool(r["selection_flip"]) for r in ledger if r["version_a"]==va and r["version_b"]==vb)})
    return ledger,summaries

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table",type=Path,default=Path("tables/gwtc_events_clean.csv"))
    p.add_argument("--v4",type=Path,default=Path("tables/gwtc_v4_frozen_structured_null.json"))
    p.add_argument("--archived-audit",type=Path,default=Path("results/gwtc_catalog_version_audit_results.json"))
    p.add_argument("--out-prefix",required=True,type=Path)
    a=p.parse_args()
    v4=json.loads(a.v4.read_text());prior=json.loads(a.archived_audit.read_text())
    if (v4.get("schema")!="GWTC_V4_STRUCTURED_NULL_FIXED_PREDICTION_V1"
        or v4.get("k_pred")!=9.7 or prior.get("schema")!="GWTC_VERSION_SENSITIVITY_V1"
        or v4.get("clean_table_sha256")!=prior.get("table_sha256")
        or v4.get("v3_frozen_sha256")!=prior.get("v3_sha256")):
        raise SystemExit("Original frozen tagged-version provenance mismatch")
    source_byte_audit=validate_frozen_csv(a.table,v4["clean_table_sha256"])
    frozen_v4_byte_audit=validate_frozen_csv(a.v4,prior["v4_sha256"])
    v=scored_versions(pd.read_csv(a.table),v4["fixed_mode"])
    ledger,summ=decompose(v)
    published={r["version"]:r for r in prior["results"]}
    for name in VERSIONS:
        measured=float(v[name]["score"].sum())
        if (not np.isclose(measured,float(published[name]["statistic_2loglr"]),rtol=0,atol=1e-8)
            or int(v[name]["selected"].sum())!=int(published[name]["positive_pastro_and_mass_n"])):
            raise SystemExit("Archived result cannot be reproduced without refit: "+name)
    oc=Path(str(a.out_prefix)+"_event_ledger.csv")
    oj=Path(str(a.out_prefix)+"_summary.json")
    if oc.exists() or oj.exists():raise SystemExit("Refusing to overwrite diagnostic")
    oc.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(ledger).to_csv(oc,index=False)
    oj.write_text(json.dumps({"schema":"GWTC_TAGGED_RELEASE_SCORE_DECOMPOSITION_V1",
       "source_table_sha256":sha(a.table),"frozen_source_csv_byte_audit":source_byte_audit,
       "frozen_v4_sha256":sha(a.v4),"frozen_v4_byte_audit":frozen_v4_byte_audit,
       "prior_version_audit_sha256":sha(a.archived_audit),
       "per_version":[{"version":ver,"selected_n":int(v[ver]["selected"].sum()),
          "score":float(v[ver]["score"].sum()),
          "reference_baseline_p_from_archived_audit":published[ver]["reference_baseline_p"]}
          for ver in VERSIONS],
       "pair_comparisons":summ,
       "scientific_status":"DESCRIPTIVE_TAGGED_VERSION_IDENTITY_DECOMPOSITION_NOT_CUMULATIVE_OR_PROSPECTIVE",
       "limitation":"Unmatched tagged entries are not physical removals; p-value differences are not additive and neither corrected for population nor independent replication."},indent=2)+"\n")
    print("Verified exact archived fixed-mode scores. Wrote",oc,oj)
    for pair in summ:
        print(pair["version_a"],"->",pair["version_b"],
              "dT=",round(pair["score_change_b_minus_a"],8),
              "shared=",round(pair["decomposition"]["shared_score_change"],8),
              "not_a=",round(pair["decomposition"]["not_listed_in_a"],8),
              "not_b=",round(pair["decomposition"]["not_listed_in_b"],8))
if __name__=="__main__":main()
