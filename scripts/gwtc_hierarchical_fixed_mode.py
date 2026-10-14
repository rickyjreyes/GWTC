#!/usr/bin/env python3
"""Fixed-mode rate-marginalized hierarchical importance-reweighting core (#8).

Requires *externally validated* whole-parameter-space population/PE-prior
and recovered-injection proposal-density ratios. It intentionally cannot
calculate them from a one-dimensional chirp-mass KDE, since that would silently
discard mass ratio, redshift, spin and detector selection. No post-hoc k scan.

For a baseline p0(theta|Lambda), define g(theta)=a cos(k ln Mchirp)+b sin(...).
At fixed Lambda and with merger rate marginalized under 1/R prior:
  log(L1/L0)=sum_i log(E_post_i[(p0/pi_PE) exp(g)]/E_post_i[p0/pi_PE])
              -N log(E_found_inj[(p0/q)exp(g)]/E_found_inj[p0/q]).
The unknown intrinsic normalization of p0 exp(g) cancels exactly. Injection
importance weights MUST already include proposal/campaign and exposure terms.
This core supplies no calibrated discovery p-value or full Bayes factor.
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
    with Path(path).open("rb") as fd:
        for chunk in iter(lambda:fd.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()
def validate_pair(ell,w,label):
    ell=np.asarray(ell,float);w=np.asarray(w,float)
    if ell.ndim!=1 or not len(ell) or ell.shape!=w.shape:
        raise ValueError(label+": missing or mismatched one-dimensional arrays")
    if not np.isfinite(ell).all() or not np.isfinite(w).all() or np.any(w<0) or not np.any(w>0):
        raise ValueError(label+": invalid log masses or importance weights")
    return ell,w
def ess(w):return float((w.sum()**2)/np.sum(w*w))
def expectation_ratio(ell,w,k,a,b):
    ell,w=validate_pair(ell,w,"importance samples")
    # g is bounded by sqrt(a^2+b^2); stable exponentiation.
    g=a*np.cos(k*ell)+b*np.sin(k*ell)
    log_r=float(logsumexp(np.log(w,where=w>0,out=np.full(w.shape,-np.inf))+g)-
                logsumexp(np.log(w,where=w>0,out=np.full(w.shape,-np.inf))))
    r=np.exp(log_r)
    # Paired delta-method estimate for log of self-normalized weighted mean.
    variance=float(np.sum(w*w*np.expm1(g-log_r)**2)/(np.sum(w)**2))
    return log_r,variance,ess(w)
def hierarchical_statistic(events,inj,k,a,b,n_generated,ess_min_event=20,ess_min_injection=None,
                           variance_limit=1.0):
    if not events or n_generated<1:raise ValueError("Need events and full generated injection count")
    il,iw=validate_pair(inj[0],inj[1],"recovered injections")
    if n_generated<len(il):raise ValueError("Total generated injections cannot be fewer than detected injections")
    xi0=float(iw.sum()/n_generated)
    if not 0<xi0<1.1:raise ValueError("Invalid baseline selection-efficiency estimate: audit proposal normalization")
    log_s,v_s,neff_s=expectation_ratio(il,iw,k,a,b)
    components=[];v_events=0
    for ell,w in events:
        logr,var,n_eff=expectation_ratio(ell,w,k,a,b)
        components.append((logr,var,n_eff))
        v_events+=var
    N=len(events)
    log_lr=float(sum(x[0] for x in components)-N*log_s)
    variance=float(v_events+N*N*v_s)
    threshold=max(20,2*N) if ess_min_injection is None else ess_min_injection
    reasons=[]
    if neff_s<threshold:reasons.append("injection_ESS_below_declared_threshold")
    if any(x[2]<ess_min_event for x in components):reasons.append("event_ESS_below_declared_threshold")
    if variance>variance_limit:reasons.append("log_likelihood_ratio_MC_variance_too_large")
    return {"statistic_2loglr":2*log_lr,"loglr_mc_se":float(np.sqrt(variance)),
      "loglr_mc_variance_estimate":variance,"event_ESS":[float(x[2]) for x in components],
      "injection_ESS":neff_s,"injection_ESS_min":threshold,"baseline_detection_fraction_estimate":xi0,
      "n_generated_injections":int(n_generated),"n_recovered_injections":len(il),
      "numerical_gate":"PASS" if not reasons else "INCOMPLETE",
      "gate_reasons":reasons,
      "scientific_status":"FIXED_HYPERPARAMETER_REWEIGHTING_CORE_NOT_POPULATION_MODEL_VALIDATION"}
def verified_npz(spec,keys):
    if not spec.get("path") or not spec.get("sha256") or sha(spec["path"]) != spec["sha256"]:
        raise ValueError("Missing/incorrect manifest file SHA256")
    with np.load(spec["path"],allow_pickle=False) as npz:
        if not set(keys).issubset(npz.files):raise ValueError("Missing required arrays in "+spec["path"])
        return [np.asarray(npz[k],float) for k in keys]
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--weights-manifest",required=True,type=Path)
    p.add_argument("--frozen",default="tables/gwtc_v4_frozen_structured_null.json",type=Path)
    p.add_argument("--output",default="tables/gwtc_hierarchical_fixed_mode_result.json",type=Path)
    args=p.parse_args()
    mf=json.loads(args.weights_manifest.read_text());fr=json.loads(args.frozen.read_text())
    if mf.get("schema")!="GWTC_HIERARCHICAL_IMPORTANCE_WEIGHTS_V1":
        raise SystemExit("Unrecognized weights manifest")
    if fr.get("schema")!="GWTC_V4_STRUCTURED_NULL_FIXED_PREDICTION_V1" or fr.get("k_pred")!=9.7:
        raise SystemExit("Incorrect frozen prediction")
    required=("population_model_and_hyperparameters","pe_reference_prior_and_reweighting",
              "injection_proposal_campaign_exposure","selection_detection_rule","official_source")
    if not all(isinstance(mf.get(key),str) and mf[key].strip() for key in required):
        raise SystemExit("Explicit population, prior, injection, detection and official-source provenance required")
    specs=mf.get("events",[])
    names=[s.get("event_name") for s in specs]
    if not names or len(names)!=len(set(names)):
        raise SystemExit("Empty/duplicate physical event manifest")
    events=[verified_npz(s,("ell","w0_over_pe_prior")) for s in specs]
    inj=verified_npz(mf["recovered_injections"],("ell_detected","w0_over_injection_proposal_detected"))
    m=fr["fixed_mode"]
    result=hierarchical_statistic(events,inj,float(m["k"]),float(m["a_cos"]),float(m["b_sin"]),
        int(mf["recovered_injections"]["n_generated"]),
        ess_min_event=int(mf.get("ess_min_event",20)),
        ess_min_injection=int(mf.get("ess_min_injection",max(20,2*len(events)))),
        variance_limit=float(mf.get("loglr_variance_limit",1.0)))
    result.update({"schema":"GWTC_HIERARCHICAL_FIXED_MODE_NUMERICAL_CORE_V1",
      "fixed_prediction_sha256":sha(args.frozen),
      "weights_manifest_sha256":sha(args.weights_manifest),
      "event_names":names,"n_physical_events":len(events),
      "evidence_interpretation":"No calibrated p or Bayes factor; validate full PE reference priors, population and injection campaign weights externally before claiming a hierarchical scientific verdict.",
      "calibrated_p":None})
    if args.output.exists():raise SystemExit("Refusing overwrite of scientific numerical result")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__=="__main__":main()
