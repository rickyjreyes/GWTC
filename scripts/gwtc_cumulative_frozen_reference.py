#!/usr/bin/env python3
"""Frozen V4/reference-KDE scores on separately supplied official cumulative tables (#12).

Only evaluate fully provenance-checked cumulative CSVs from #12's manifest;
NEVER infer cumulative snapshots from incremental GWOSC tagged rows. These are
conditional reference-KDE p-values, not population/injection-corrected,
independent, or prospective replication tests. No frozen parameter refits.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from gwtc_catalog_version_audit import VERSIONS
from gwtc_cumulative_manifest_audit import load_releases
from gwtc_source_provenance import validate_frozen_csv
from gwtc_unbinned_kde import fixed_log_likelihood_ratio,fixed_model_null

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def evaluate_cumulative(tables,v3,v4,null_n,seed):
    if null_n<1000:raise ValueError("Require >=1000 predeclared reference nulls")
    m=v4["fixed_mode"];rows=[]
    for i,version in enumerate(VERSIONS):
        df=tables[version]
        p=pd.to_numeric(df.p_astro,errors="coerce")
        z=pd.to_numeric(df.M_chirp,errors="coerce")
        cut=p>=.5
        good=cut & np.isfinite(z) & (z>0)
        ell=np.log(z.loc[good].to_numpy(float))
        row={"version":version,"canonical_primary_rows":int(len(df)),
             "pastro_pass_rows":int(cut.sum()),"positive_pastro_and_mass_n":int(good.sum()),
             "missing_or_nonpositive_mass_after_pastro_cut":int((cut&~good).sum()),
             "frozen_k":float(m["k"]),"null_n":int(null_n),"seed":int(seed+i*1000003),
             "reference_calibration":"SINGLE_FIXED_V4_MODE_VS_FROZEN_V3_TRAINING_KDE",
             "population_or_selection_corrected":False,"prospective":False}
        if len(ell):
            observed=float(fixed_log_likelihood_ratio(ell,m["k"],m["a_cos"],m["b_sin"],m["normalizer_Z"]))
            null=fixed_model_null(kde_centers=v3["train_ell"],bandwidth=v3["selected_bandwidth"],
                 n_events=len(ell),k=m["k"],a=m["a_cos"],b=m["b_sin"],
                 normalizer=m["normalizer_Z"],null_n=null_n,seed=seed+i*1000003)
            exceed=int((null>=observed).sum())
            row.update({"statistic_2loglr":observed,"null_exceed":exceed,
                 "fixed_reference_p":float((1+exceed)/(1+null_n))})
        else:
            row.update({"statistic_2loglr":None,"null_exceed":None,"fixed_reference_p":None})
        rows.append(row)
    return rows

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest",required=True,type=Path)
    p.add_argument("--v3",type=Path,default=Path("tables/gwtc_v3_frozen_unbinned_kde_mode.json"))
    p.add_argument("--v4",type=Path,default=Path("tables/gwtc_v4_frozen_structured_null.json"))
    p.add_argument("--archived-audit",type=Path,default=Path("results/gwtc_catalog_version_audit_results.json"))
    p.add_argument("--null-n",type=int,default=10000)
    p.add_argument("--seed",type=int,default=20261014)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    if a.output.exists():raise SystemExit("Refusing to overwrite cumulative result")
    manifest=json.loads(a.manifest.read_text());tables=load_releases(manifest)
    v3=json.loads(a.v3.read_text());v4=json.loads(a.v4.read_text())
    prior=json.loads(a.archived_audit.read_text())
    if (v4.get("schema")!="GWTC_V4_STRUCTURED_NULL_FIXED_PREDICTION_V1" or v4.get("k_pred")!=9.7
        or prior.get("v3_sha256")!=v4.get("v3_frozen_sha256")
        or prior.get("schema")!="GWTC_VERSION_SENSITIVITY_V1"):
        raise SystemExit("Frozen reference provenance mismatch")
    v3_provenance=validate_frozen_csv(a.v3,v4["v3_frozen_sha256"])
    v4_provenance=validate_frozen_csv(a.v4,prior["v4_sha256"])
    rows=evaluate_cumulative(tables,v3,v4,a.null_n,a.seed)
    result={"schema":"GWTC_OFFICIAL_CUMULATIVE_FIXED_REFERENCE_V1",
      "source_manifest_sha256":sha(a.manifest),
      "frozen_v3_byte_audit":v3_provenance,"frozen_v4_byte_audit":v4_provenance,
      "sources":[{"version":s["version"],"source_url":s["source_url"],
                  "table_sha256":s["sha256"],"published_utc":s["published_utc"]}
                 for s in manifest["releases"]],
      "results":rows,"scientific_status":"CONDITIONAL_CUMULATIVE_REFERENCE_SENSITIVITY_NOT_INDEPENDENT_OR_HIERARCHICAL",
      "limitation":"Verified table hashes and stated source URLs do not independently certify official cumulative completeness or waveform equivalence. Overlapping cumulative catalogs are dependent; reference p-values exclude population, PE and injection selection uncertainty."}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+"\n")
    print("Wrote frozen cumulative reference comparison:",a.output)
    for r in rows:
        print(r["version"],"n=",r["positive_pastro_and_mass_n"],
              "score=",r["statistic_2loglr"],"fixed_KDE_reference_p=",r["fixed_reference_p"])
if __name__=="__main__":main()
