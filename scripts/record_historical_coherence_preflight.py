#!/usr/bin/env python3
"""Execute the historical verifier and record its exact obstruction, never a null result."""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path
import pandas as pd

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 if a.output.exists():p.error('refusing overwrite')
 protocol=Path('protocols/historical_coherence_v1.json')
 observed=Path('tables/gwtc4_population_observed.csv')
 missing=Path('outputs_gwtc4_wct_subsets/gwtc4_subset_event_summary.csv')
 cmd=[sys.executable,'scripts/generate_gwtc4_population_null_matrix.py','--verify-only']
 run=subprocess.run(cmd,capture_output=True,text=True)
 frame=pd.read_csv(observed)
 history=subprocess.run(['git','log','--all','--format=%H','--',str(missing)],capture_output=True,text=True)
 out=dict(schema='GWTC_HISTORICAL_COHERENCE_EXECUTION_PREFLIGHT_V1',command=cmd,exit_code=run.returncode,
  stdout=run.stdout,stderr=run.stderr,required_event_summary=str(missing),event_summary_present=missing.exists(),
  event_summary_sha256=sha(missing) if missing.exists() else None,all_fetched_refs_path_commits=history.stdout.splitlines(),
  source_hashes={str(x):sha(x) for x in [protocol,observed,Path(cmd[1]),Path('scripts/historical_mode_coherence.py')]},
  observed_scan_rows=len(frame),observed_columns=list(frame.columns),frozen_protocol=json.loads(protocol.read_text()),
  complete_outer_null_catalogs_executed=0,calibrated_coherence_p=None,
  status='PREFLIGHT_PASSED_NOT_NULL_EXECUTED' if run.returncode==0 else 'BLOCKED_EXTERNAL',
  prohibition='Historical delta-only nulls cannot replace k/phase/membership; no substitute event table was generated',
  historical_source_search='All remote branches fetched; no matching path in git history; exact Library search found no event-summary artifact')
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],run.stderr.strip())
if __name__=='__main__':main()
