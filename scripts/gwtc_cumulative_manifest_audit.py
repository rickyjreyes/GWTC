#!/usr/bin/env python3
"""Official *cumulative* GWTC release reconciliation (#12); no prospective claim.

Distinct from the old tagged-incremental-row sensitivity diagnostic. Every
cumulative table must be supplied by the operator with official source URL,
SHA256, official release ID/date and canonical physical event identity. An
absent event is NOT silently counted as a physically removed merger. Rename
mappings are explicit, externally sourced declarations, never fuzzy matched.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from gwtc_catalog_version_audit import FIELDS

VERSIONS=("GWTC-4.0","GWTC-4.1","GWTC-5.0")
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load_releases(manifest):
    if manifest.get("schema")!="GWTC_CUMULATIVE_MANIFEST_V1":
        raise ValueError("Wrong official cumulative manifest schema")
    specs=manifest.get("releases",[])
    if tuple(x.get("version") for x in specs)!=VERSIONS:
        raise ValueError("Require three exact ordered release IDs")
    tables={}
    for spec in specs:
        if not all(isinstance(spec.get(k),str) and spec[k].strip() for k in
                   ("source_url","published_utc","path","sha256")):
            raise ValueError("Missing official release provenance")
        if not spec["source_url"].startswith("https://"):raise ValueError("Invalid source URL")
        if sha(spec["path"])!=spec["sha256"]:raise ValueError("Cumulative table SHA mismatch")
        df=pd.read_csv(spec["path"])
        needed={"commonName","is_primary_entry","M_chirp","p_astro"}|set(FIELDS)
        if not needed.issubset(df):raise ValueError("Incomplete cumulative table fields")
        if df.commonName.isna().any():raise ValueError("Null merger identity")
        primary=df.is_primary_entry.astype(str).str.strip().str.lower().isin(("true","1","yes"))
        df=df.loc[primary].copy()
        if df.commonName.duplicated().any():raise ValueError("Ambiguous canonical identities")
        tables[spec["version"]]=df.set_index("commonName",drop=False)
    return tables

def reconcile(tables, aliases=()):
    """Explicit renames only; no ungrounded removal/new-event interpretation."""
    if not aliases:aliases=[]
    seen=set()
    for a in aliases:
        if not all(isinstance(a.get(k),str) and a[k].strip() for k in
                   ("version_a","version_b","from_name","to_name","official_identity_evidence_url")):
            raise ValueError("Every rename requires complete external official identity evidence")
        key=(a["version_a"],a["version_b"],a["from_name"])
        if key in seen:raise ValueError("Duplicate rename mapping")
        seen.add(key)
    rows=[];summary=[]
    for va,vb in combinations(VERSIONS,2):
        a=tables[va];b=tables[vb]
        renames={x["from_name"]:x["to_name"] for x in aliases if x["version_a"]==va and x["version_b"]==vb}
        if len(set(renames.values()))!=len(renames):raise ValueError("Many-to-one alias requires manual reconciliation")
        if any(k not in a.index or v not in b.index for k,v in renames.items()):
            raise ValueError("Rename maps absent official primary identity")
        if any(k in b.index or v in a.index for k,v in renames.items()):
            raise ValueError("Ambiguous rename conflicts with unchanged primary identity")
        common=set(a.index)&set(b.index)
        linked=[(v,v) for v in sorted(common)]+sorted(renames.items())
        claimed_a={x for x,_ in linked};claimed_b={y for _,y in linked}
        for x,y in linked:
            changes=[]
            for f in FIELDS:
                left=pd.to_numeric(pd.Series([a.loc[x,f]]),errors="coerce").iloc[0]
                right=pd.to_numeric(pd.Series([b.loc[y,f]]),errors="coerce").iloc[0]
                if pd.isna(left) and pd.isna(right):continue
                if pd.isna(left)!=pd.isna(right) or not np.isclose(left,right,rtol=0,atol=1e-12):
                    changes.append(f)
            rows.append({"version_a":va,"version_b":vb,"name_a":x,"name_b":y,
                "status":"OFFICIAL_RENAME_DECLARED" if x!=y else ("REVISED_FIELDS" if changes else "UNCHANGED_FIELDS"),
                "changed_fields":";".join(changes)})
        for x in sorted(set(a.index)-claimed_a):
            rows.append({"version_a":va,"version_b":vb,"name_a":x,"name_b":"",
                 "status":"NOT_LISTED_IN_B_NOT_PHYSICAL_REMOVAL","changed_fields":""})
        for y in sorted(set(b.index)-claimed_b):
            rows.append({"version_a":va,"version_b":vb,"name_a":"","name_b":y,
                 "status":"NOT_LISTED_IN_A_NEW_TO_THIS_SNAPSHOT_NOT_PROSPECTIVE","changed_fields":""})
        pair=[r for r in rows if r["version_a"]==va and r["version_b"]==vb]
        summary.append({"version_a":va,"version_b":vb,"canonical_primary_a":len(a),
            "canonical_primary_b":len(b),
            "unmatched_a":sum(r["status"].startswith("NOT_LISTED_IN_B") for r in pair),
            "unmatched_b":sum(r["status"].startswith("NOT_LISTED_IN_A") for r in pair),
            "declared_renames":len(renames),
            "revised_matched":sum(bool(r["changed_fields"]) for r in pair)})
    return rows,summary

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest",required=True,type=Path)
    p.add_argument("--out-prefix",required=True,type=Path)
    args=p.parse_args()
    mf=json.loads(args.manifest.read_text())
    tables=load_releases(mf)
    rows,summ=reconcile(tables,mf.get("declared_renames",[]))
    out_csv=Path(str(args.out_prefix)+"_event_diff.csv")
    out_json=Path(str(args.out_prefix)+"_summary.json")
    if out_csv.exists() or out_json.exists():raise SystemExit("Refuse overwrite")
    out_csv.parent.mkdir(parents=True,exist_ok=True)
    pd.DataFrame(rows).to_csv(out_csv,index=False)
    out_json.write_text(json.dumps({"schema":"GWTC_CUMULATIVE_RECONCILIATION_V1",
      "manifest_sha256":sha(args.manifest),
      "inputs":[{"version":x["version"],"source_url":x["source_url"],"sha256":x["sha256"],"published_utc":x["published_utc"]} for x in mf["releases"]],
      "comparisons":summ,"scientific_status":"DESCRIPTIVE_OFFICIAL_IDENTITY_AUDIT_NOT_PREDICTIVE_OR_POPULATION_CALIBRATED",
      "reference_KDE_tagged_row_diagnostic":"results/gwtc_catalog_version_audit_results.json",
      "missing_identity_interpretation":"An event absent from a catalog is not proof that a physical merger was removed."},indent=2)+"\n")
    print("Wrote",out_csv,out_json)
if __name__=="__main__":main()
