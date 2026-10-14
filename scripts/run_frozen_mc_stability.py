#!/usr/bin/env python3
"""Re-run frozen V1/V3/V4 nulls with a PREDECLARED seed/depth matrix (#14).

Writes new outputs to a fresh directory. Never overwrites original artifacts
or freezes. Reports exact/binomial Clopper-Pearson intervals for exceedance
fraction, plus add-one Monte Carlo p and between-run classification stability.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import beta
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def clopper_pearson(k,n,confidence=.95):
    if not 0<=k<=n or n<1:raise ValueError("Invalid exceedance count")
    alpha=1-confidence
    low=0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1))
    high=1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k))
    return low,high
def extract_rows(mode,csv):
    df=pd.read_csv(csv)
    if mode=="V1":
        return [("fixed_polynomial",float(df.iloc[0].holdout_delta_deviance),
                 int(df.iloc[0].null_n),float(df.iloc[0].holdout_p))]
    if mode=="V3":
        return [("fixed_KDE",float(df.iloc[0].delta_2logl_holdout),
                 int(df.iloc[0].null_n),float(df.iloc[0].p_holdout))]
    return [(str(row.scenario),float(row.delta_2logl_holdout),
             int(row.null_n),float(row.p_structured_null)) for row in df.itertuples(index=False)]
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan",type=Path,default=ROOT/"protocols/mc_stability_v1.json")
    p.add_argument("--out-dir",type=Path,required=True)
    p.add_argument("--modes",nargs="+",default=["V1","V3","V4"],choices=("V1","V3","V4"))
    args=p.parse_args()
    plan=json.loads(args.plan.read_text())
    if plan["schema"]!="GWTC_FROZEN_MC_STABILITY_V1":raise SystemExit("Wrong frozen stability plan")
    if args.out_dir.exists():raise SystemExit("Choose a NEW output directory; refuse overwrite")
    args.out_dir.mkdir(parents=True)
    table=ROOT/"tables/gwtc_events_clean.csv"
    manifest=ROOT/"tables/gwtc_holdout_manifest.csv"
    before={"plan":sha(args.plan),"table":sha(table),"manifest":sha(manifest)}
    rows=[];orig_hash={}
    for mode in args.modes:
        spec=plan["modes"][mode];frozen=ROOT/spec["frozen_path"]
        orig_hash[mode]=sha(frozen)
        for seed in plan["seeds"][mode]:
            for n in plan["depths"]:
                target=args.out_dir/f"{mode}_seed{seed}_n{n}.csv"
                cmd=[sys.executable,str(ROOT/spec["script"]),"--table",str(table),
                     "--manifest",str(manifest),spec["frozen_arg"],str(frozen),
                     "--seed",str(seed),"--null-n",str(n),"--output",str(target)]
                subprocess.run(cmd,cwd=ROOT,check=True)
                if sha(frozen)!=orig_hash[mode]:raise RuntimeError("Frozen model changed DURING stability run")
                for scenario,observed,actual_n,pval in extract_rows(mode,target):
                    if actual_n!=n:raise RuntimeError("Evaluator changed requested null depth")
                    k=int(round(pval*(n+1)-1))
                    if abs(pval-(1+k)/(n+1))>1e-8:raise RuntimeError("Invalid add-one empirical p")
                    lo,hi=clopper_pearson(k,n,plan["confidence"])
                    rows.append({"mode":mode,"scenario":scenario,"seed":seed,"null_n":n,
                        "observed_statistic":observed,"exceedances":k,"empirical_p":pval,
                        "binomial_interval_low":lo,"binomial_interval_high":hi,
                        "class_pass_point":bool(observed>0 and pval<plan["alpha"]),
                        "source_result_sha256":sha(target)})
    if before!={"plan":sha(args.plan),"table":sha(table),"manifest":sha(manifest)}:
        raise RuntimeError("Provenance changed during numerical-stability runs")
    df=pd.DataFrame(rows);df.to_csv(args.out_dir/"matrix.csv",index=False)
    status=[]
    for (mode,scenario),sub in df.groupby(["mode","scenario"],sort=True):
        if not np.allclose(sub.observed_statistic,sub.observed_statistic.iloc[0],rtol=0,atol=1e-7):
            verdict="MODEL_OR_DATA_DRIFT"
        elif sub.class_pass_point.nunique()!=1:verdict="SEED_OR_DEPTH_CLASSIFICATION_UNSTABLE"
        elif sub.class_pass_point.iloc[0] and (sub.binomial_interval_high<plan["alpha"]).all():
            verdict="NUMERICALLY_STABLE_PASS"
        elif (not sub.class_pass_point.iloc[0]) and (sub.binomial_interval_low>=plan["alpha"]).all():
            verdict="NUMERICALLY_STABLE_NONPASS"
        else:verdict="MC_PRECISION_UNRESOLVED"
        status.append({"mode":mode,"scenario":scenario,"stability_verdict":verdict,
             "p_min":float(sub.empirical_p.min()),"p_max":float(sub.empirical_p.max()),
             "exceedance_counts":sub.exceedances.tolist()})
    (args.out_dir/"summary.json").write_text(json.dumps({
       "schema":"GWTC_MC_STABILITY_RESULTS_V1","frozen_hashes":orig_hash,
       "input_hashes":before,"plan_sha256":before["plan"],"status":status,
       "interpretation":"Seed/depth sensitivity of identical frozen models; never choose best seed; no new independent data."
    },indent=2)+"\n")
    print("Saved independent matrix and status at",args.out_dir)
if __name__=="__main__":main()
