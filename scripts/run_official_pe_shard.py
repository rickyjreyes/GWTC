#!/usr/bin/env python3
"""Bounded raw-file retrieval and lossless source-chirp projections for CI.

Each original file is size/MD5 checked and SHA256 hashed before extracting all
samples. Raw files are removed only after their receipt and verified projection
are saved. These projections are local derivatives, not official HDF5 products.
"""
import argparse,csv,hashlib,json,subprocess,time
from pathlib import Path
import h5py,numpy as np
from recover_official_pe_inputs import verify_file
from gwtc_posterior_uncertainty import file_sha

def project(row,raw,directory):
 receipt=verify_file(row,raw)
 with h5py.File(raw) as f:mass=np.asarray(f[row['dataset_path']].fields('chirp_mass_source')[:],dtype='<f8')
 out=directory/(row['event_name']+'.hdf5');values=np.empty(len(mass),dtype=[('chirp_mass_source','<f8')]);values['chirp_mass_source']=mass
 with h5py.File(out,'x') as f:
  f.create_dataset('posterior_samples',data=values,compression='gzip')
  f.attrs['derived_product']='lossless projection of all official source-chirp posterior samples'
  f.attrs['raw_sha256']=receipt['sha256'];f.attrs['official_md5']=row['official_md5'];f.attrs['waveform']=row['waveform']
  f.attrs['unit']='M_sun';f.attrs['frame']='source';f.attrs['weighting']='equal posterior draws'
 with h5py.File(out) as f:
  if not np.array_equal(mass,f['posterior_samples']['chirp_mass_source']):raise ValueError('projection mismatch')
 receipt.update(projection_filename=out.name,projection_sha256=file_sha(out),
  source_mass_array_sha256=hashlib.sha256(mass.tobytes()).hexdigest(),projection_lossless_verified=True)
 return receipt

def accepted_projection(row,receipt,directory):
 """Only reuse a projection bound to this exact official identity and bytes."""
 for key in ('event_name','source_url','official_md5','bytes','dataset_path','waveform'):
  if receipt.get(key)!=row.get(key):raise ValueError('cached source identity mismatch: '+key)
 if receipt.get('status')!='VERIFIED' or not receipt.get('md5_verified') or not receipt.get('projection_lossless_verified'):
  raise ValueError('cached raw verification missing')
 path=directory/receipt['projection_filename']
 if path.name!=row['event_name']+'.hdf5' or path.parent!=directory:raise ValueError('unexpected projection filename')
 if file_sha(path)!=receipt['projection_sha256']:raise ValueError('cached projection checksum mismatch')
 return receipt

def retrieve(row,raw):
 """Resume partial transfers, with bounded retries for TLS resets/timeouts."""
 url=row['source_url'].replace('/api/records/','/records/').removesuffix('/content')+'?download=1'
 errors=[]
 for attempt in range(4):
  p=subprocess.run(['curl','--http1.1','-L','--fail','--silent','--show-error',
   '--connect-timeout','60','--max-time','900','--speed-time','120','--speed-limit','512',
   '--continue-at','-',url,'-o',str(raw)],capture_output=True,text=True)
  if not p.returncode:return
  errors.append(p.stderr.strip())
  if p.returncode==33:raw.unlink(missing_ok=True) # server explicitly refused Range
  if attempt<3:time.sleep(5*(attempt+1))
 raise RuntimeError('\n'.join(errors))

def run(mapping,shard,shards,output,seed=None):
 output.mkdir(parents=True,exist_ok=True);scratch=output/'raw';scratch.mkdir(exist_ok=True)
 obj=json.loads(mapping.read_text());allrows=sorted(obj['events'],key=lambda x:x['event_name'])
 rows=allrows[shard::shards];receipts=[]
 for row in rows:
  raw=scratch/row.get('filename',row['event_name'])
  try:
   if row['status']!='MAPPED_NOT_DOWNLOADED':raise ValueError('metadata mapping not accepted')
   cached=seed/(row['event_name']+'.receipt.json') if seed else None
   if cached and cached.exists():
    r=accepted_projection(row,json.loads(cached.read_text()),seed)
    (output/r['projection_filename']).write_bytes((seed/r['projection_filename']).read_bytes())
   else:
    retrieve(row,raw);r=project(row,raw,output)
   receipts.append(r)
   (output/(row['event_name']+'.receipt.json')).write_text(json.dumps(r,indent=2)+'\n')
   raw.unlink(missing_ok=True)
  except (OSError,ValueError,KeyError,RuntimeError) as exc:
   receipts.append({**row,'status':'FAILED','error':str(exc)})
   raw.unlink(missing_ok=True)
  print(row['event_name'],receipts[-1]['status'],flush=True)
  (output/f'shard{shard}.json').write_text(json.dumps({'mapping_sha256':file_sha(mapping),'events':receipts},indent=2)+'\n')
 if any(r['status']!='VERIFIED' for r in receipts):raise SystemExit(2)

def collect(directory,mapping):
 receipts=[json.loads(p.read_text()) for p in sorted(directory.glob('*.receipt.json'))]
 expected=json.loads(mapping.read_text())['events'];identities={r['event_name']:r for r in expected}
 byname={r['event_name']:accepted_projection(identities[r['event_name']],r,directory) for r in receipts}
 if len(byname)!=len(receipts):raise ValueError('duplicate shard event')
 manifest=[]
 for name,r in sorted(byname.items()):
  path=directory/r['projection_filename']
  if file_sha(path)!=r['projection_sha256'] or not r['projection_lossless_verified'] or not r['md5_verified']:raise ValueError('projection or raw source verification failed')
  manifest.append(dict(event_name=name,path=str(path),dataset_path='posterior_samples',sha256=r['projection_sha256'],source_url=r['source_url'],weight_path=''))
 with (directory/'posterior_manifest.csv').open('x') as f:
  w=csv.DictWriter(f,fieldnames=['event_name','path','dataset_path','sha256','source_url','weight_path']);w.writeheader();w.writerows(manifest)
 summary={'schema':'GWTC_OFFICIAL_PE_PROJECTION_INVENTORY_V1','mapping_sha256':file_sha(mapping),'expected':len(expected),'verified':len(byname),
  'missing':sorted({r['event_name'] for r in expected}-set(byname)),'raw_files_verified_then_discarded':True,'events':receipts}
 (directory/'source_receipts.json').write_text(json.dumps(summary,indent=2)+'\n')
 print('OFFICIAL_SOURCE_COVERAGE',len(byname),'/',len(expected),flush=True)
 if summary['missing']:raise SystemExit(2)

def main():
 p=argparse.ArgumentParser();p.add_argument('--mapping',type=Path,default=Path('results/official_2026-10-15/pe_input_mapping.json'))
 p.add_argument('--shard',type=int,default=0);p.add_argument('--shards',type=int,default=8);p.add_argument('--output',type=Path,required=True);p.add_argument('--seed',type=Path);p.add_argument('--collect',action='store_true');a=p.parse_args()
 if a.collect:collect(a.output,a.mapping)
 else:run(a.mapping,a.shard,a.shards,a.output,a.seed)
if __name__=='__main__':main()
