#!/usr/bin/env python3
"""Read-only Nyquist annotation for archived binned exploratory GWTC tables (#15).

Uses the historical identity n_star = k_best * active log width / (2*pi);
hence k_N=pi*B/width. Does NOT recompute p, change legacy verdict, or
overwrite the frozen historical table. Inferred width is labelled explicitly.
Unbinned results and historical 168 continuous-x scans are outside this gate.
"""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd

def annotate_legacy(frame):
    required={"bin_count","k_best","n_star"}
    if not required.issubset(frame):raise ValueError("Missing binned legacy columns")
    out=frame.copy()
    b=pd.to_numeric(out.bin_count,errors="coerce").to_numpy(float)
    k=pd.to_numeric(out.k_best,errors="coerce").to_numpy(float)
    n=pd.to_numeric(out.n_star,errors="coerce").to_numpy(float)
    ok=np.isfinite(b)&np.isfinite(k)&np.isfinite(n)&(b>=2)&(k>0)&(n>0)
    width=np.full(len(out),np.nan);width[ok]=2*np.pi*n[ok]/k[ok]
    kn=np.full(len(out),np.nan);kn[ok]=np.pi*b[ok]/width[ok]
    status=np.full(len(out),"UNDETERMINED_MISSING_LEGACY_FIELDS",dtype=object)
    status[ok&(k>=kn*(1-1e-10))]="ALIASED_AT_OR_ABOVE_NYQUIST"
    status[ok&(k<kn*(1-1e-10))]="RESOLVED_BELOW_NYQUIST"
    out["legacy_implied_active_log_width"]=width
    out["nyquist_k"]=kn
    out["nyquist_status"]=status
    out["original_verdict_preserved"]=out["verdict_label"] if "verdict_label" in out else "NOT_RECORDED"
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    args=p.parse_args()
    if args.input.resolve()==args.output.resolve() or args.output.exists():
        raise SystemExit("Refusing to overwrite historical artifact")
    result=annotate_legacy(pd.read_csv(args.input))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    result.to_csv(args.output,index=False)
    print("Historical table SHA256:",hashlib.sha256(args.input.read_bytes()).hexdigest())
    print(result.nyquist_status.value_counts().to_string())
    print("Read-only annotations written:",args.output)
if __name__=="__main__":main()
