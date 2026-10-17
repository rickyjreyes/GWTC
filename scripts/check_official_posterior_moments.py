#!/usr/bin/env python3
"""Exact empirical-posterior moments for the additive frozen statistic.

This checks the posterior-draw Monte Carlo without estimating a null p-value.
Weights are normalized within each verified event; event independence is the
same measurement-diagnostic assumption used by the frozen sampler.
"""
import argparse,json,math
from pathlib import Path
import h5py,numpy as np
from run_official_pe_shard import accepted_projection
from gwtc_posterior_uncertainty import file_sha

def main():
 p=argparse.ArgumentParser();p.add_argument('--receipts',type=Path,required=True);p.add_argument('--diagnostic',type=Path,required=True)
 p.add_argument('--mapping',type=Path,default=Path('results/official_2026-10-15/pe_input_mapping.json'))
 p.add_argument('--frozen',type=Path,default=Path('tables/gwtc_v4_frozen_structured_null.json'));p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 if a.output.exists():p.error('refusing overwrite')
 d=json.loads(a.diagnostic.read_text());m=json.loads(a.frozen.read_text())['fixed_mode'];source=json.loads(a.mapping.read_text())['events']
 if d['coverage']!='COMPLETE' or d['frozen_sha256']!=file_sha(a.frozen) or len(source)!=d['input_events']:
  raise ValueError('incomplete or changed frozen diagnostic')
 means=[];variances=[]
 for row in sorted(source,key=lambda x:x['event_name']):
  r=accepted_projection(row,json.loads((a.receipts/(row['event_name']+'.receipt.json')).read_text()),a.receipts)
  with h5py.File(a.receipts/r['projection_filename']) as f:ell=np.log(f['posterior_samples']['chirp_mass_source'])
  t=2*(m['a_cos']*np.cos(m['k']*ell)+m['b_sin']*np.sin(m['k']*ell)-np.log(m['normalizer_Z']))
  means.append(float(t.mean()));variances.append(float(t.var()))
 mean=math.fsum(means);var=math.fsum(variances);se=math.sqrt(var/d['draws']);delta=d['posterior_draw_score_mean']-mean
 if not math.isclose(mean,math.fsum(r['posterior_mean_score_contribution'] for r in d['events']),rel_tol=0,abs_tol=1e-10):
  raise ValueError('posterior event ledger disagrees with direct empirical moments')
 out=dict(schema='GWTC_EMPIRICAL_POSTERIOR_MOMENT_CHECK_V1',diagnostic_sha256=file_sha(a.diagnostic),mapping_sha256=file_sha(a.mapping),
  events=len(source),empirical_posterior_mean_statistic=mean,empirical_posterior_variance=var,
  empirical_posterior_standard_deviation=math.sqrt(var),draws=d['draws'],draw_mean=d['posterior_draw_score_mean'],
  monte_carlo_standard_error_of_mean=se,draw_mean_minus_exact_mean=delta,mc_standard_error_units=delta/se,
  point_statistic=d['point_statistic_on_matched_events'],relative_exact_mean_attenuation=1-mean/d['point_statistic_on_matched_events'],
  calibrated_p=None,interpretation='Empirical posterior expectation of additive log statistic, not log of an integrated likelihood or a hierarchical Bayes factor')
 a.output.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
