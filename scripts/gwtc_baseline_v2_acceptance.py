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

ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
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
