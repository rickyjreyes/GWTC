#!/usr/bin/env python3
"""Conditional baseline-hyperparameter mixture of fixed-mode hierarchical LR (#8).

For weights P(Lambda | d,H0) fixed INDEPENDENTLY of the periodic model,
B_10 = sum_Lambda P(Lambda|d,H0) * L1(d|Lambda)/L0(d|Lambda).
This does NOT itself fit p0, validate event PE/reference priors, infer
injection proposals, marginalize rate under a different rate prior, or
calibrate a null p-value. Inputs are outputs of gwtc_hierarchical_fixed_mode.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""):h.update(block)
    return h.hexdigest()

def conditional_mixture(items):
    if not items or len({i["model_id"] for i in items})!=len(items):
        raise ValueError("Empty/duplicate population models")
    weights=np.asarray([i["baseline_posterior_weight"] for i in items],float)
    stats=np.asarray([i["statistic_2loglr"] for i in items],float)
    if not (np.isfinite(weights).all() and np.isfinite(stats).all() and
            (weights>0).all() and np.isclose(weights.sum(),1,rtol=0,atol=1e-8)):
        raise ValueError("Unnormalized/invalid baseline posterior weights or LR")
    if any(i.get("numerical_gate")!="PASS" for i in items):
        raise ValueError("Numerically incomplete model cannot enter mixture")
    if len({tuple(i["event_names"]) for i in items})!=1 or len({i["fixed_prediction_sha256"] for i in items})!=1:
        raise ValueError("Mismatched event set or frozen mode across models")
    logbf=float(logsumexp(np.log(weights)+stats/2))
    return {"conditional_model_averaged_2loglr":2*logbf,
       "individual_2loglr":stats.tolist(),
       "model_ids":[i["model_id"] for i in items],
       "baseline_posterior_weights":weights.tolist(),
       "fixed_prediction_sha256":items[0]["fixed_prediction_sha256"],
       "event_names":items[0]["event_names"],"calibrated_p":None,
       "scientific_status":"CONDITIONAL_ON_EXTERNALLY_VALIDATED_BASELINE_POSTERIOR_AND_INJECTION_WEIGHTS_NOT_A_CALIBRATED_BAYES_FACTOR"}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--grid-manifest",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    mf=json.loads(args.grid_manifest.read_text())
    if mf.get("schema")!="GWTC_HIERARCHICAL_BASELINE_GRID_V1" or mf.get("weights_frozen_under_no_mode_model_before_test") is not True:
        raise SystemExit("Require frozen H0 baseline posterior mixture provenance")
    if not all(isinstance(mf.get(x),str) and mf[x].strip() for x in ("posterior_weight_derivation","population_prior_provenance","injection_selection_consistency")):
        raise SystemExit("Missing full scientific provenance")
    items=[]
    for spec in mf["models"]:
        path=Path(spec["result_path"])
        if not spec.get("sha256") or sha(path)!=spec["sha256"]:raise SystemExit("Input result SHA mismatch")
        row=json.loads(path.read_text())
        if row.get("schema")!="GWTC_HIERARCHICAL_FIXED_MODE_NUMERICAL_CORE_V1":
            raise SystemExit("Incorrect numerical result schema")
        row["model_id"]=spec["model_id"]
        row["baseline_posterior_weight"]=float(spec["baseline_posterior_weight"])
        items.append(row)
    out=conditional_mixture(items)
    out.update({"schema":"GWTC_HIERARCHICAL_CONDITIONAL_MIXTURE_V1",
       "grid_manifest_sha256":sha(args.grid_manifest),
       "individual_result_sha256":[m["sha256"] for m in mf["models"]]})
    if args.output.exists():raise SystemExit("Refusing to overwrite comparison")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out,indent=2))
if __name__=="__main__":main()
