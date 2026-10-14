#!/usr/bin/env python3
"""Historical mode recurrence across overlapping subset scans (#10).

Scoring rule must be committed BEFORE any complete-catalog outer-null details
are generated. The existing null-Delta/p matrices alone are insufficient:
they do not contain peak k, phase, or dynamic subset event membership.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

REQUIRED=("scan_id","x_var","y_var","events","k_best","phase","n_obs","delta_x")
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def checked(frame):
    if not set(REQUIRED).issubset(frame.columns):
        raise ValueError("Missing full coherence k/phase/n/membership fields")
    if frame.scan_id.duplicated().any():raise ValueError("Duplicate scan_id")
    for key in ("k_best","phase","n_obs","delta_x"):
        if not np.isfinite(pd.to_numeric(frame[key],errors="coerce")).all():
            raise ValueError("Nonfinite peak parameter "+key)
    if (frame.delta_x<=0).any():raise ValueError("Invalid active x span")
    ncalc=frame.k_best*frame.delta_x/(2*math.pi)
    if not np.allclose(frame.n_obs,ncalc,atol=1e-7,rtol=1e-7):
        raise ValueError("Inconsistent n_obs vs k and active domain")
    return frame

def coherence_statistic(frame, config):
    frame=checked(frame)
    unique=[]; seen=set()
    for r in frame.sort_values("scan_id").itertuples(index=False):
        ev=frozenset(e for e in str(r.events).split(",") if e)
        if not ev:raise ValueError("Missing physical-event membership")
        identity=(str(r.x_var),str(r.y_var),ev)
        if identity in seen:continue # Duplicate labels are NOT replications
        seen.add(identity)
        unique.append((r,ev))
    num=0.;den=0.;pairs=0
    for i,(a,ae) in enumerate(unique):
        for b,be in unique[i+1:]:
            if (a.x_var,a.y_var)!=(b.x_var,b.y_var) or ae==be:continue
            jaccard=len(ae&be)/len(ae|be)
            if jaccard<config["min_jaccard"]:continue
            phase_diff=abs(math.atan2(math.sin(a.phase-b.phase),math.cos(a.phase-b.phase)))
            coherent=(abs(a.k_best-b.k_best)<=config["k_tolerance"] and
                      phase_diff<=config["phase_tolerance_radians"])
            den+=jaccard
            num+=jaccard*float(coherent)
            pairs+=1
    if not pairs or den<=0:raise ValueError("No eligible distinct, overlapping subset pairs")
    return {"score":float(num/den),"pair_count":pairs,"deduplicated_scans":len(unique)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--observed",type=Path,default=Path("tables/gwtc4_population_observed.csv"))
    p.add_argument("--null-details",required=True,type=Path)
    p.add_argument("--config",type=Path,default=Path("protocols/historical_coherence_v1.json"))
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    cfg=json.loads(args.config.read_text())
    if cfg["schema"]!="GWTC4_HISTORICAL_COHERENCE_FROZEN_V1":raise SystemExit("Wrong frozen protocol")
    obs=pd.read_csv(args.observed)
    if len(obs)!=168:raise SystemExit("Historical observed requires all 168 original rows")
    if "scan_id" not in obs.columns:
        obs["scan_id"]=[f"scan_{i:03d}" for i in range(len(obs))]
    checked(obs)
    obs_score=coherence_statistic(obs,cfg)
    null=pd.read_csv(args.null_details)
    if "replicate" not in null:raise SystemExit("Null detail matrix lacks replicate IDs")
    scores=[]
    for rid,group in null.groupby("replicate",sort=True):
        if set(group.scan_id)!=set(obs.scan_id) or len(group)!=len(obs):
            raise SystemExit(f"Incomplete full-catalog null replicate {rid}")
        same=group.set_index("scan_id").loc[obs.scan_id]
        if not ((same.x_var.to_numpy()==obs.x_var.to_numpy()).all()
                and (same.y_var.to_numpy()==obs.y_var.to_numpy()).all()):
            raise SystemExit("Null changed frozen scan-pair identity")
        scores.append(coherence_statistic(group,cfg)["score"])
    if len(scores)<2:raise SystemExit("Insufficient outer-null catalogs")
    n_ge=int(np.sum(np.asarray(scores)>=obs_score["score"]))
    pval=float((1+n_ge)/(1+len(scores)))
    result={"schema":"GWTC4_HISTORICAL_COHERENCE_RESULT_V1",
      "observed_sha256":sha(args.observed),"null_details_sha256":sha(args.null_details),
      "config_sha256":sha(args.config),"observed":obs_score,
      "outer_catalogs":len(scores),"outer_exceed":n_ge,"global_p":pval,
      "verdict":"PARTIAL_OUTER_MC_DEPTH" if len(scores)<1000 else
                ("GLOBALLY_RARE_UNDER_THIS_DECLARED_NULL" if pval<.01 else "NOT_GLOBALLY_RARE"),
      "interpretation":"Exploratory historical recurrence under stated complete-catalog outer-null; not an independent frozen GWTC future prediction."}
    if args.output.exists():raise SystemExit("Refusing result overwrite")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
