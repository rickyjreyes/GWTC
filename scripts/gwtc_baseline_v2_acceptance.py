#!/usr/bin/env python3
"""Original baseline-v2 (#3) execution gate, preserving PARTIAL/Class III.

Machine-readable acceptance audit. A posterior-draw diagnostic is necessary,
but never sufficient for a Class-I promotion: historical leakage/provenance,
hierarchical population/selection, and unseen prospective replication remain
independent limitations. Does not rewrite V1/V3/V4 historical results.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import math
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def archived_consistency(root=ROOT):
    """Check archived cross-file identities; never reinterpret the original tests."""
    files={name:root/"tables"/filename for name,filename in {
        "cv":"gwtc_baseline_v2_cv.csv","v2":"gwtc_frozen_mode.json",
        "v3":"gwtc_v3_frozen_unbinned_kde_mode.json",
        "v2_result":"gwtc_frozen_holdout_result.csv",
        "v3_result":"gwtc_v3_unbinned_kde_holdout_result.csv",
        "manifest":"gwtc_holdout_manifest.csv","table":"gwtc_events_clean.csv"}.items()}
    absent=sorted(k for k,v in files.items() if not v.exists())
    if absent:return {"status":"INCOMPLETE_ARCHIVED_INPUTS","missing":absent}
    cv=pd.read_csv(files["cv"]);v2=json.loads(files["v2"].read_text())
    v3=json.loads(files["v3"].read_text())
    a=pd.read_csv(files["v2_result"]).iloc[0]
    b=pd.read_csv(files["v3_result"]).iloc[0]
    held=pd.read_csv(files["manifest"])
    selected=cv.loc[cv.selected.astype(str).str.lower().eq("true") &
                    cv.fold.astype(str).eq("mean")]
    if len(selected)!=1 or selected.iloc[0].data_split!="train" or selected.iloc[0].baseline_family!="polynomial_poisson" or int(selected.iloc[0].degree)!=int(v2["baseline_degree"]):
        raise ValueError("Archived V2 CV selected row/frozen degree mismatch")
    if v3.get("bandwidth_selection")!="training-only leave-one-out log likelihood" or v3.get("baseline_family")!="gaussian_kde":
        raise ValueError("Missing training-selected non-polynomial baseline")
    if v3.get("manifest_sha256")!=sha(files["manifest"]) or v3.get("clean_table_sha256")!=sha(files["table"]):
        raise ValueError("Frozen V3 input hashes do not match archived artifacts")
    if not math.isclose(float(v2["k_star"]),float(a["k_frozen"]),rel_tol=0,abs_tol=1e-10) or not math.isclose(float(v3["k_star"]),float(b["k_frozen"]),rel_tol=0,abs_tol=1e-10):
        raise ValueError("Archived fixed-mode evaluation has been retuned")
    if int(a.training_n)!=173 or int(a.holdout_positive_n)!=104 or int(b.holdout_n_positive_finite)!=104:
        raise ValueError("Frozen source population/count drift")
    names=held.loc[held.split=="holdout","event_name"].astype(str)
    if names.duplicated().any() or len(names)!=int(b.holdout_manifest_rows):
        raise ValueError("Invalid archived holdout manifest")
    if str(b.scientific_status)!="ROBUSTNESS_TEST_NOT_INDEPENDENT_REPLICATION" or str(a.class_i_claim).lower() not in ("false","0"):
        raise ValueError("Historical promotion/robustness status drift")
    width=float(v2["bin_edges"][-1])-float(v2["bin_edges"][0])
    k_n=math.pi*int(v2["bins"])/width
    if float(v2["k_star"])>=k_n*(1-1e-10) or not math.isclose(k_n,float(v2["nyquist_k"]),rel_tol=1e-10):
        raise ValueError("Archived V2 mode fails sampling Nyquist")
    return {"status":"ARCHIVED_CONSISTENT_NOT_PROSPECTIVE_REPLICATION",
       "training_n":173,"holdout_positive_n":104,"v2_selected_polynomial_degree":int(v2["baseline_degree"]),
       "v3_baseline":"training-only Gaussian KDE",
       "v2_frozen_k":float(v2["k_star"]),"v3_frozen_k":float(v3["k_star"]),
       "v2_reference_p":float(a.holdout_p),"v3_reference_p":float(b.p_holdout),
       "v2_nyquist_k":k_n,
       "interpretation":"V3 is a robustness reanalysis of previously inspected GWTC-5. The earlier pre-strict CV had holdout exposure. No new independent replication."}

def baseline_gate(posterior=None,posterior_sha=None,root=ROOT):
    frozen=root/"tables/gwtc_v4_frozen_structured_null.json"
    required={
      "training_cv":root/"tables/gwtc_baseline_v2_cv.csv",
      "training_selected_nonpolynomial_KDE":root/"tables/gwtc_v3_frozen_unbinned_kde_mode.json",
      "unbinned_holdout":root/"tables/gwtc_v3_unbinned_kde_holdout_result.csv",
      "frozen_holdout_manifest":root/"tables/gwtc_holdout_manifest.csv",
      "fixed_mode_no_rescan_holdout":root/"tables/gwtc_frozen_holdout_result.csv",
      "historical_negative_outcome_registry":root/"tables/gwtc4_population_global_result.json",
    }
    gates={k:{"status":"ARCHIVED_REQUIRES_PROVENANCE_REVIEW" if f.is_file() else "MISSING",
              "path":str(f.relative_to(root)),"sha256":sha(f) if f.is_file() else None}
           for k,f in required.items()}
    post={"status":"NOT_RUN_NO_OFFICIAL_POSTERIOR_DIAGNOSTIC"}
    if posterior is not None:
        pp=Path(posterior)
        if not posterior_sha or sha(pp)!=posterior_sha:raise ValueError("Posterior result SHA mismatch")
        obj=json.loads(pp.read_text())
        if obj.get("schema")!="GWTC_V4_POSTERIOR_SENSITIVITY_V1" or obj.get("frozen_sha256")!=sha(frozen):
            raise ValueError("Wrong frozen V4 posterior sensitivity result")
        ok=(obj.get("coverage")=="COMPLETE" and obj.get("input_events")==obj.get("holdout_expected_positive_mass")
            and not obj.get("missing_event_names") and obj.get("calibrated_p") is None)
        post={"status":"EXECUTED_COMPLETE_COVERAGE_DIAGNOSTIC_REQUIRES_SCIENTIFIC_REVIEW" if ok else
             "INCOMPLETE_MEASUREMENT_UNCERTAINTY_GATE",
             "result_sha256":posterior_sha,"n_events":obj.get("input_events"),
             "n_expected":obj.get("holdout_expected_positive_mass"),
             "mean_minus_point":obj.get("posterior_mean_minus_point_statistic"),
             "posterior_draw_quantiles":obj.get("posterior_draw_score_quantiles_05_50_95")}
    gates["posterior_measurement_uncertainty"]=post
    gates["archived_cross_file_consistency"]=archived_consistency(root)
    return {"schema":"GWTC_BASELINE_V2_ACCEPTANCE_RECONCILIATION_V1",
       "gates":gates,"original_holdout_provenance_caveat":"Earlier pre-strict CV had already seen holdout; do not claim an untouched independent historical replication.",
       "classification":"PARTIAL_CLASS_III_NOT_PROMOTED",
       "promotion":False,"unresolved_dependencies":["#8 full selection/population challenge","#18 unseen prospective replication"],
       "scientific_verdict":"ORIGINAL_ACCEPTANCE_NOT_CLOSED" if post["status"]!="EXECUTED_COMPLETE_COVERAGE_DIAGNOSTIC_REQUIRES_SCIENTIFIC_REVIEW" else
          "POSTERIOR_DIAGNOSTIC_DONE_REQUIRE_SCIENTIFIC_INTERPRETATION_AND_ORIGINAL_GATE_AUDIT"}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--posterior-result",type=Path)
    p.add_argument("--posterior-sha256")
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    if a.output.exists():raise SystemExit("Refusing to overwrite acceptance record")
    out=baseline_gate(a.posterior_result,a.posterior_sha256)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,indent=2)+"\n")
    print(out["scientific_verdict"],"=>",a.output)
if __name__=="__main__":main()
