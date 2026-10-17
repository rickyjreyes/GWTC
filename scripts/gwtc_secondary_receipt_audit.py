#!/usr/bin/env python3
"""Expose official raw-field availability from checksum-bound ingestion receipts.

The saved source-chirp projections deliberately do not contain the secondary
values. Therefore this inventories the *original* verified dataset field names;
it does not certify secondary missingness, units, producer version or analysis.
"""
import argparse,json
from pathlib import Path
import pandas as pd
from run_official_pe_shard import accepted_projection
from audit_official_gwtc_observables import audit_field_names
from gwtc_posterior_uncertainty import file_sha

def audit(mapping,directory):
 rows=[]
 for source in mapping['events']:
  path=directory/(source['event_name']+'.receipt.json')
  receipt=accepted_projection(source,json.loads(path.read_text()),directory) if path.exists() else None
  fields=audit_field_names(receipt['posterior_field_names']) if receipt else {}
  for observable in ('E_rad','peak_luminosity'):
   f=fields.get(observable,{'status':'RAW_SOURCE_NOT_RETRIEVED','exact_fields_found':[]})
   rows.append(dict(event_name=source['event_name'],observable=observable,status=f['status'],
    field_names=';'.join(f['exact_fields_found']),source_url=source['source_url'],declared_dataset=source['dataset_path'],
    file_sha256=receipt['sha256'] if receipt else '',receipt_sha256=file_sha(path) if receipt else '',
    official_md5=source['official_md5'],raw_checksum_verified=bool(receipt),
    units_verified=False,frame_convention_verified=False,secondary_values_retained=False,analysis_ready=False))
 return pd.DataFrame(rows)

def main():
 p=argparse.ArgumentParser();p.add_argument('--mapping',type=Path,default=Path('results/official_2026-10-15/pe_input_mapping.json'))
 p.add_argument('--receipts',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 if a.output.exists():p.error('refusing overwrite')
 rows=audit(json.loads(a.mapping.read_text()),a.receipts);rows.to_csv(a.output,index=False)
 print(rows.groupby(['observable','status']).size().to_string())
if __name__=='__main__':main()
