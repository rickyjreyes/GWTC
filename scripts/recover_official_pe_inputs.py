#!/usr/bin/env python3
"""Resolve official preferred-waveform identities before PE ingestion (#11/#13).

The prepare step reads only catalog/PE metadata. Download verifies the entire
official file size and MD5, computes SHA256, and inventories its true fields.
Neither incomplete bytes nor missing units count as accepted posterior data.
"""
import argparse,concurrent.futures,hashlib,json,re,subprocess
from pathlib import Path
from urllib.parse import urlparse
import h5py,numpy as np,pandas as pd
from gwtc_posterior_uncertainty import file_sha
from gwtc_source_provenance import validate_frozen_csv
from audit_official_gwtc_observables import audit_field_names

def preferred_event(obj,name):
 events=obj.get('events',{})
 found=[v for v in events.values() if v.get('commonName')==name]
 if len(found)!=1:raise ValueError('ambiguous official event identity')
 event=found[0]
 pe=[v for v in event.get('parameters',{}).values() if v.get('pipeline_type')=='pe' and v.get('is_preferred') is True]
 if len(pe)!=1:raise ValueError('ambiguous preferred official PE analysis')
 p=pe[0];url=p.get('data_url','');parsed=urlparse(url)
 if parsed.scheme!='https' or parsed.hostname!='zenodo.org':raise ValueError('unsupported official PE source host')
 if not re.search(r'/'+re.escape(name)+r'[-/]',url.replace('-'+name+'-','/'+name+'-')):raise ValueError('event/file identity mismatch')
 if p.get('chirp_mass_source_unit')!='M_sun' or not np.isfinite(p.get('chirp_mass_source',np.nan)):
  raise ValueError('unverified source-frame solar-mass units')
 return dict(event_name=name,waveform=p['waveform_family'],source_url=url,
  dataset_path=p['waveform_family']+'/posterior_samples',mass_field='chirp_mass_source',
  unit='M_sun',frame='source',official_point_mass=p['chirp_mass_source'],
  official_date_added=p.get('date_added'),catalog_version=event['catalog.shortName'])

def fetch_json(url,path):
 if not path.exists():
  r=subprocess.run(['curl','-L','--fail','--silent','--show-error','--connect-timeout','15','--max-time','60','--retry','1',url,'-o',str(path)],capture_output=True,text=True)
  if r.returncode:
   path.unlink(missing_ok=True);raise RuntimeError(r.stderr.strip())
 return json.loads(path.read_text())

def prepare(directory,records,output):
 directory.mkdir(parents=True,exist_ok=True)
 freeze=json.loads(Path('tables/gwtc_v4_frozen_structured_null.json').read_text())
 validate_frozen_csv('tables/gwtc_events_clean.csv',freeze['clean_table_sha256'])
 validate_frozen_csv('tables/gwtc_holdout_manifest.csv',freeze['manifest_sha256'])
 table=pd.read_csv('tables/gwtc_events_clean.csv');h=pd.read_csv('tables/gwtc_holdout_manifest.csv')
 held=set(h.loc[h.split=='holdout','event_name']);primary=table.is_primary_entry.astype(str).str.lower().isin(['true','1','yes'])
 selected=table[primary&table.commonName.isin(held)&np.isfinite(table.M_chirp)&(table.M_chirp>0)]
 catalog=json.loads((records/'GWTC-5.0.json').read_text())['events']
 cat={v['commonName']:v for v in catalog.values()}
 products={}
 for path in sorted(records.glob('zenodo*_record.json')):
  record=json.loads(path.read_text())
  for name,entry in record['files']['entries'].items():
   products[entry['links']['content']]={'filename':name,'bytes':entry['size'],'official_md5':entry['checksum'].removeprefix('md5:'),
     'record_id':record['id'],'record_version_index':record['versions']['index'],'record_snapshot_sha256':file_sha(path)}
 def one(row):
  name=row.commonName;url=cat[name]['jsonurl'];path=directory/(name+'.json')
  try:
   obj=fetch_json(url,path);d=preferred_event(obj,name)
   if d['source_url'] not in products:raise ValueError('preferred PE file absent from pinned official release inventory')
   if float(d['official_point_mass'])!=float(row.M_chirp):raise ValueError('official preferred point mass changed from frozen table')
   return {**d,**products[d['source_url']],'metadata_url':url,'metadata_sha256':file_sha(path),'status':'MAPPED_NOT_DOWNLOADED'}
  except (OSError,ValueError,RuntimeError,KeyError) as exc:return {'event_name':name,'metadata_url':url,'status':'BLOCKED','error':str(exc)}
 with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(one,selected.sort_values('commonName').itertuples()))
 out={'schema':'GWTC_OFFICIAL_PE_RECOVERY_V1','expected_events':len(selected),'events':rows,
  'units_policy':'official preferred PE source chirp mass in M_sun, exactly paired with frozen point estimate',
  'no_posterior_statistic_executed':True,'records_sha256':{p.name:file_sha(p) for p in sorted(records.glob('zenodo*_record.json'))}}
 with output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
 print('mapped',sum(r['status']=='MAPPED_NOT_DOWNLOADED' for r in rows),'/',len(rows),flush=True)

def verify_file(row,path):
 if path.stat().st_size!=row['bytes']:raise ValueError('full official byte size mismatch')
 md5=hashlib.md5();sha=hashlib.sha256()
 with path.open('rb') as f:
  for data in iter(lambda:f.read(8*1024*1024),b''):md5.update(data);sha.update(data)
 if md5.hexdigest()!=row['official_md5']:raise ValueError('official MD5 mismatch')
 with h5py.File(path) as f:
  dataset=f[row['dataset_path']];fields=dataset.dtype.names
  if not fields or row['mass_field'] not in fields:raise ValueError('source-frame posterior field missing')
  masses=np.asarray(dataset.fields(row['mass_field'])[:],float)
  if len(masses)<2 or not np.isfinite(masses).all() or np.any(masses<=0):raise ValueError('invalid source-frame posterior samples')
  weighting_fields=[x for x in fields if x.lower() in ['weight','weights','log_weight','log_weights']]
  if weighting_fields:raise ValueError('explicit posterior weights require independent normalization policy')
  receipt={**row,'path':str(path),'sha256':sha.hexdigest(),'md5_verified':True,
   'posterior_n':len(masses),'weight_ess':len(masses),'weight_policy':'equal normalized posterior draws; no pooling across events',
   'posterior_field_names':list(fields),'secondary_fields':audit_field_names(fields),
   'secondary_analysis_ready':False,'status':'VERIFIED'}
 return receipt

def download(manifest,directory,output,event_names=None):
 directory.mkdir(parents=True,exist_ok=True);m=json.loads(manifest.read_text());rows=[]
 for row in m['events']:
  if event_names and row['event_name'] not in event_names:continue
  if row['status']!='MAPPED_NOT_DOWNLOADED':rows.append(row);continue
  path=directory/row['filename']
  try:
   if not path.exists():
    url=row['source_url'].replace('/api/records/','/records/').removesuffix('/content')+'?download=1'
    partial=path.with_suffix('.partial')
    p=subprocess.run(['curl','-L','--fail','--silent','--show-error','--connect-timeout','20','--max-time','1800','-C','-',url,'-o',str(partial)],capture_output=True,text=True)
    if p.returncode:raise RuntimeError(p.stderr.strip())
    partial.rename(path)
   rows.append(verify_file(row,path))
  except (OSError,ValueError,RuntimeError,KeyError) as exc:rows.append({**row,'status':'BLOCKED_EXTERNAL_TRANSFER_OR_METADATA','error':str(exc)})
  print(row['event_name'],rows[-1]['status'],flush=True)
 with output.open('x') as f:json.dump({'schema':'GWTC_PE_DOWNLOAD_RECEIPTS_V1','events':rows},f,indent=2);f.write('\n')

def main():
 p=argparse.ArgumentParser();p.add_argument('--prepare',action='store_true');p.add_argument('--directory',type=Path,required=True)
 p.add_argument('--records',type=Path,default=Path('results/official_2026-10-15'));p.add_argument('--manifest',type=Path)
 p.add_argument('--output',type=Path,required=True);p.add_argument('--event',action='append');a=p.parse_args()
 if a.output.exists():p.error('refusing overwrite')
 if a.prepare:prepare(a.directory,a.records,a.output)
 else:download(a.manifest,a.directory,a.output,a.event)
if __name__=='__main__':main()
