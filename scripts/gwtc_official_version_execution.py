#!/usr/bin/env python3
"""Execute the existing frozen comparator on pinned official incremental rows.

No union is represented as an official historical cumulative snapshot. GPS
collisions/rename candidates are exposed, never silently merged or deleted.
"""
import argparse,concurrent.futures,json
from pathlib import Path
from itertools import combinations
import numpy as np,pandas as pd
from gwtc_catalog_version_audit import VERSIONS,FIELDS,version_diff
from gwtc_cumulative_frozen_reference import evaluate_cumulative
from gwtc_version_score_decomposition import scored_versions,decompose
from gwtc_source_provenance import validate_frozen_csv
from recover_official_pe_inputs import fetch_json
from gwtc_posterior_uncertainty import file_sha

def normalize(obj,version):
 rows=[];names=set();gps={}
 for key,v in obj['events'].items():
  name=v['commonName'];g=float(v['GPS'])
  if name in names:raise ValueError('duplicate official identity')
  if v['catalog.shortName']!=version or key!=name+'-v'+str(v['version']):raise ValueError('catalog/event version mismatch')
  if not np.isfinite(g):raise ValueError('invalid GPS')
  if g in gps:raise ValueError('ambiguous same-GPS physical identities')
  names.add(name);gps[g]=name
  mass=v.get('chirp_mass_source')
  if mass is not None and (not np.isfinite(mass) or mass<=0 or v.get('chirp_mass_source_unit')!='M_sun'):
   raise ValueError('invalid official source chirp mass or units')
  row=dict(commonName=name,catalog_version=version,official_entry=key,event_version=v['version'],GPS=g,
   metadata_url=v['jsonurl'],M_chirp=mass,M_total=v.get('total_mass_source'),M_final=v.get('final_mass_source'),
   chirp_mass_field='chirp_mass_source',chirp_mass_unit=v.get('chirp_mass_source_unit'),mass_frame='source')
  row.update({k:v.get(k) for k in FIELDS if k not in row});rows.append(row)
 return pd.DataFrame(rows).sort_values('commonName').reset_index(drop=True)

def reconcile(a,b):
 shared=set(a.commonName)&set(b.commonName)
 ai=a.set_index('commonName');bi=b.set_index('commonName')
 for n in shared:
  if ai.loc[n,'GPS']!=bi.loc[n,'GPS']:raise ValueError('shared event has changed physical GPS identity')
 aa=a[~a.commonName.isin(shared)];bb=b[~b.commonName.isin(shared)]
 possible=aa[['commonName','GPS']].merge(bb[['commonName','GPS']],on='GPS',suffixes=('_a','_b')).to_dict('records')
 return dict(shared=len(shared),not_listed_in_a=len(bb),not_listed_in_b=len(aa),
  same_gps_rename_candidates=possible,physical_removal_inferred=False,automatic_alias_merges=0)

def pe_record(row,directory):
 path=directory/(row.commonName+'.json')
 obj=fetch_json(row.metadata_url,path)
 candidates=[e for e in obj['events'].values() if e.get('commonName')==row.commonName]
 if len(candidates)!=1:raise ValueError('ambiguous metadata identity')
 event=candidates[0]
 if event['catalog.shortName']!=row.catalog_version:raise ValueError('wrong metadata catalog')
 pe=[p for p in event.get('parameters',{}).values() if p.get('pipeline_type')=='pe' and p.get('is_preferred') is True]
 if len(pe)!=1:raise ValueError('preferred PE is ambiguous/missing')
 p=pe[0]
 if p.get('chirp_mass_source_unit')!='M_sun' or float(p['chirp_mass_source'])!=float(row.M_chirp):
  raise ValueError('preferred PE differs from catalog source-mass definition/value')
 return dict(version=row.catalog_version,event_name=row.commonName,metadata_url=row.metadata_url,
  metadata_sha256=file_sha(path),preferred_waveform=p['waveform_family'],source_url=p.get('data_url'),
  date_added=p.get('date_added'),mass_field='chirp_mass_source',units='M_sun',frame='source',status='VERIFIED_METADATA')

def main():
 p=argparse.ArgumentParser();p.add_argument('--sources',type=Path,default=Path('results/official_2026-10-15'))
 p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
 tables={v:normalize(json.loads((a.sources/(v+'.json')).read_text()),v) for v in VERSIONS}
 pe=[]
 for v,table in tables.items():
  directory=a.sources/('event_metadata_'+v);directory.mkdir(exist_ok=True)
  selected=table[(table.p_astro>=.5)&table.M_chirp.notna()]
  def one(row):
   try:return pe_record(row,directory)
   except (ValueError,OSError,KeyError,RuntimeError) as exc:return dict(version=v,event_name=row.commonName,status='BLOCKED',error=str(exc))
  with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:pe.extend(pool.map(one,selected.itertuples()))
  table.to_csv(a.output/(v+'_manifest.csv'),index=False)
  print(v,'PE verified',sum(x['status']=='VERIFIED_METADATA' and x['version']==v for x in pe),'/',len(selected),flush=True)
 v3p=Path('tables/gwtc_v3_frozen_unbinned_kde_mode.json');v4p=Path('tables/gwtc_v4_frozen_structured_null.json')
 v3=json.loads(v3p.read_text());v4=json.loads(v4p.read_text());prior=json.loads(Path('results/gwtc_catalog_version_audit_results.json').read_text())
 provenance={'v3':validate_frozen_csv(v3p,v4['v3_frozen_sha256']),'v4':validate_frozen_csv(v4p,prior['v4_sha256'])}
 if v4['k_pred']!=9.7:raise ValueError('wrong frozen mode')
 scores=evaluate_cumulative(tables,v3,v4,10000,20261014)
 for r in scores:r['release_scope']='OFFICIAL_INCREMENTAL_TAGGED_ROWS_NOT_CUMULATIVE'
 diffs=[];identity=[]
 for va,vb in combinations(VERSIONS,2):
  diffs+=version_diff(tables[va].set_index('commonName',drop=False),tables[vb].set_index('commonName',drop=False),va,vb)
  identity.append(dict(version_a=va,version_b=vb,**reconcile(tables[va],tables[vb])))
 pd.DataFrame(diffs).to_csv(a.output/'event_diff.csv',index=False)
 ledger,ds=decompose(scored_versions(pd.concat(tables.values()),v4['fixed_mode']))
 pd.DataFrame(ledger).to_csv(a.output/'score_change_ledger.csv',index=False)
 out=dict(schema='GWTC_OFFICIAL_INCREMENTAL_VERSION_EXECUTION_V1',sources={v:{'url':'https://gwosc.org/eventapi/json/'+v+'/',
  'sha256':file_sha(a.sources/(v+'.json'))} for v in VERSIONS},frozen_provenance=provenance,
  results=scores,identity_reconciliation=identity,additive_score_decomposition=ds,preferred_PE_metadata=pe,
  cumulative_snapshot_status='BLOCKED: versioned authoritative all-run as-of snapshots not supplied; current incremental rows cannot substitute',
  scientific_status='RETROSPECTIVE_CATALOG_SENSITIVITY_NOT_PROSPECTIVE_OR_POPULATION_CALIBRATED',
  limitations=['PE source mass conventions verified from preferred metadata; waveform families can differ',
   'Reference-KDE null ignores population, event PE uncertainty and detection selection',
   'GWTC-4.0 portal default PE was updated after initial release; this is an as-retrieved snapshot, not original publication bytes'])
 (a.output/'results.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(scores,indent=2),flush=True)
 if any(x['status']!='VERIFIED_METADATA' for x in pe):raise SystemExit(2)
if __name__=='__main__':main()
