#!/usr/bin/env python3
"""Fail-closed official-field definition preflight for secondary observables (#13).

Consumes the #13 exact-field coverage audit and a separately reviewed official
definition manifest. It DOES NOT infer units from names, compute E_rad from
mass differences, run unregistered scans, or claim any physical measurement.
Only matched, operator-attested definitions are conditionally eligible for a
new frozen exploratory-family protocol; all other cases remain blocked.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def preflight(coverage, definitions):
    if definitions.get("schema")!="GWTC_SECONDARY_OFFICIAL_DEFINITIONS_V1":
        raise ValueError("Wrong official definition manifest")
    if coverage.duplicated(["event_name","observable"]).any():
        raise ValueError("Duplicate event-observable source rows")
    entries=definitions.get("observables",[])
    if len({x.get("observable") for x in entries})!=len(entries):
        raise ValueError("Duplicate official definition records")
    registered={x.get("observable"):x for x in entries}
    allowed={"E_rad","peak_luminosity"}
    if not set(registered).issubset(allowed):raise ValueError("Unregistered secondary variable")
    rows=[]
    for row in coverage.itertuples(index=False):
        d=registered.get(row.observable,{})
        reason=[]
        if row.status!="PRESENT_EXACT_NAME":reason.append("OFFICIAL_EXACT_FIELD_ABSENT_OR_AMBIGUOUS")
        if d.get("exact_field") not in str(row.field_names).split(";"):
            reason.append("FIELD_NOT_IN_SOURCE_DATASET")
        for k in ("official_definition_url","release_id","definition_sha256","units","frame_convention","exact_field"):
            if not isinstance(d.get(k),str) or not d[k].strip():reason.append("MISSING_"+k.upper())
        if d.get("independent_definition_reviewed") is not True:
            reason.append("OFFICIAL_DEFINITION_REVIEW_NOT_ATTESTED")
        if d.get("official_definition_url") and not d["official_definition_url"].startswith("https://"):
            reason.append("INVALID_SOURCE_URL")
        rows.append({"event_name":row.event_name,"observable":row.observable,
            "file_sha256":row.file_sha256,"field_names":row.field_names,
            "preflight_status":"CONDITIONAL_DEFINITION_REVIEWED_NOT_SCIENTIFIC_RESULT" if not reason else "BLOCKED",
            "reasons":";".join(sorted(set(reason))),"units":d.get("units",""),
            "frame_convention":d.get("frame_convention",""),
            "definition_url":d.get("official_definition_url","")})
    return pd.DataFrame(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--coverage",required=True,type=Path)
    p.add_argument("--definitions",required=True,type=Path)
    p.add_argument("--out-prefix",required=True,type=Path)
    a=p.parse_args()
    defs=json.loads(a.definitions.read_text())
    result=preflight(pd.read_csv(a.coverage,keep_default_na=False),defs)
    outcsv=Path(str(a.out_prefix)+"_rows.csv")
    outjson=Path(str(a.out_prefix)+"_summary.json")
    if outcsv.exists() or outjson.exists():raise SystemExit("Refuse overwrite")
    outcsv.parent.mkdir(parents=True,exist_ok=True)
    result.to_csv(outcsv,index=False)
    outjson.write_text(json.dumps({"schema":"GWTC_SECONDARY_PREFLIGHT_V1",
       "coverage_sha256":sha(a.coverage),"definitions_sha256":sha(a.definitions),
       "blocked_rows":int((result.preflight_status=="BLOCKED").sum()),
       "conditional_rows":int((result.preflight_status!="BLOCKED").sum()),
       "scientific_status":"NO_SCANS_OR_NULL_CALIBRATION_PERFORMED",
       "next_gate":"Predeclare transforms, missingness, version compatibility and complete exploratory null family before analysis.",
       "calibrated_p":None},indent=2)+"\n")
    print("Wrote fail-closed secondary preflight:",outcsv,outjson)
if __name__=="__main__":main()
