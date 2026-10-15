#!/usr/bin/env python3
"""Separate explanatory figure for the ARCHIVED binned exploratory Nyquist audit.

Never silently mixes the aliased B20 peak into frozen V2/V3 evidence plots.
Consumes read-only overlay CSVs; does not rewrite old tables or publication art.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
DEFAULT=ROOT/"results"/"legacy_nyquist"/"gwtc_bin_stress_annotated.csv"

def validate_overlay(df):
    required={"bin_count","k_best","nyquist_k","nyquist_status"}
    if not required.issubset(df) or df.bin_count.duplicated().any() or df.empty:
        raise ValueError("Incomplete or duplicate legacy bin annotation")
    for r in df.itertuples(index=False):
        if not (r.bin_count>=2 and r.k_best>0 and r.nyquist_k>0):
            raise ValueError("Invalid legacy peak/Nyquist")
        expected=("ALIASED_AT_OR_ABOVE_NYQUIST" if r.k_best>=r.nyquist_k*(1-1e-10)
                  else "RESOLVED_BELOW_NYQUIST")
        if r.nyquist_status!=expected:raise ValueError("Contradictory saved Nyquist classification")
    return df.sort_values("bin_count").copy()

def make_figure(df,output):
    df=validate_overlay(df)
    out=Path(output)
    if out.exists():raise FileExistsError("Never overwrite publication audit figure")
    if out.suffix.lower() not in (".svg",".png",".pdf"):raise ValueError("Use SVG, PNG or PDF")
    fig,ax=plt.subplots(figsize=(9,5))
    ax.plot(df.bin_count,df.nyquist_k,marker="o",label="Per-binning Nyquist limit")
    first_alias=True;first_below=True
    for r in df.itertuples(index=False):
        alias=r.nyquist_status.startswith("ALIASED")
        label=("Aliased historical peak" if first_alias else None) if alias else ("Below-Nyquist historical peak" if first_below else None)
        ax.scatter([r.bin_count],[r.k_best],marker="x" if alias else "o",s=90,label=label)
        if alias: first_alias=False
        else: first_below=False
        ax.annotate("ALIASED: not a resolved mode" if alias else "Below Nyquist",
                    (r.bin_count,r.k_best),xytext=(5,7),textcoords="offset points",fontsize=8)
    ax.set_xlabel("Archived exploratory bin count B")
    ax.set_ylabel("Angular log-frequency k")
    ax.set_title("Read-only historical bin-stress sampling audit (not primary evidence)")
    ax.text(.02,-.24,"Original p-values/verdicts remain archival, not alias-corrected inference. "
            "Frozen V2 B=40 and V3 unbinned tests are separate.",
            transform=ax.transAxes,fontsize=8,va="top")
    ax.legend(frameon=False,fontsize=9)
    fig.tight_layout()
    out.parent.mkdir(parents=True,exist_ok=True)
    if out.suffix.lower()==".svg":fig.savefig(out,metadata={"Date":None},bbox_inches="tight")
    else:fig.savefig(out,bbox_inches="tight")
    plt.close(fig)
    return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--overlay",type=Path,default=DEFAULT)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    print("Wrote separate archival sampling figure:",make_figure(pd.read_csv(a.overlay),a.output))
if __name__=="__main__":main()
