#!/usr/bin/env python3
"""Plot a complete official posterior diagnostic; never turn its spread into a p-value."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
 p=argparse.ArgumentParser();p.add_argument('--diagnostic',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 d=json.loads(a.diagnostic.read_text());rows=d['events']
 if d['coverage']!='COMPLETE' or d['missing_event_names'] or len(rows)!=d['holdout_expected_positive_mass']:
  raise SystemExit('Refuse full-holdout figure for incomplete coverage')
 point=np.array([r['point_score_contribution'] for r in rows]);mean=np.array([r['posterior_mean_score_contribution'] for r in rows])
 if not np.isfinite(point).all() or not np.isfinite(mean).all() or not np.isclose(point.sum(),d['point_statistic_on_matched_events'],rtol=0,atol=1e-10):
  raise SystemExit('Inconsistent per-event score ledger')
 if a.output.exists():raise SystemExit('Refusing figure overwrite')
 q=np.array(d['posterior_draw_score_quantiles_05_50_95'])
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.hashsalt':'GWTC-official-posterior-v1'})
 fig,(ax,bx)=plt.subplots(1,2,figsize=(11,4.8),gridspec_kw={'width_ratios':[1.1,1]})
 ax.axvline(0,color='#b6bcc3',lw=1)
 ax.scatter([d['point_statistic_on_matched_events']],[1],color='#285a8e',s=75,zorder=3,label='Catalog point estimate')
 ax.errorbar([q[1]],[0],xerr=[[q[1]-q[0]],[q[2]-q[1]]],fmt='o',color='#b04e32',capsize=5,lw=3,label='Draw median and 5–95% interval')
 ax.scatter([d['posterior_draw_score_mean']],[0],marker='D',s=42,color='#252525',zorder=4,label='Draw mean')
 ax.set(yticks=[0,1],yticklabels=['Posterior draws','Point values'],ylim=(-.7,1.55),xlabel='Frozen statistic T = 2 log LR',title='One posterior draw per event')
 ax.legend(loc='upper left',bbox_to_anchor=(0,-.2),frameon=False,fontsize=8)
 lo=min(point.min(),mean.min());hi=max(point.max(),mean.max());pad=.07*(hi-lo)
 bx.plot([lo-pad,hi+pad],[lo-pad,hi+pad],color='#aeb4bc',lw=1,ls='--',label='Equal contribution')
 bx.scatter(point,mean,s=24,alpha=.65,color='#285a8e',edgecolor='none')
 bx.set(xlim=(lo-pad,hi+pad),ylim=(lo-pad,hi+pad),xlabel='Point-estimate event contribution',ylabel='Posterior-mean event contribution',title=f'Paired contributions: {len(rows)} events')
 fig.suptitle('Official GWTC-5.0 measurement-uncertainty diagnostic',fontweight='bold',y=.98)
 fig.text(.5,.015,'Frozen k and phase • posterior spread, not a calibrated p-value • no hierarchical population/selection correction',ha='center',fontsize=9,color='#525a65')
 fig.subplots_adjust(left=.12,right=.98,top=.85,bottom=.28,wspace=.35)
 a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,metadata={'Date':None},bbox_inches='tight');plt.close(fig)
 print(a.output)
if __name__=='__main__':main()
