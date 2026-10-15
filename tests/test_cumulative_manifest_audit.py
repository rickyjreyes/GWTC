import sys
from pathlib import Path
import hashlib
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from gwtc_cumulative_manifest_audit import load_releases,reconcile
from gwtc_catalog_version_audit import FIELDS

def mock_row(name,m):
    row={"commonName":name,"is_primary_entry":True,"p_astro":.9}
    for f in FIELDS:row[f]=float("nan")
    row["M_chirp"]=m
    return row
def test_distinct_cumulative_manifest_verifies_hash_and_never_infers_physical_removal(tmp_path):
    rows=[[mock_row("a",10),mock_row("b",20)],
          [mock_row("a",11),mock_row("renamed",20)],
          [mock_row("a",11),mock_row("renamed",20),mock_row("new",30)]]
    specs=[]
    for version,rr in zip(("GWTC-4.0","GWTC-4.1","GWTC-5.0"),rows):
        p=tmp_path/(version+".csv");pd.DataFrame(rr).to_csv(p,index=False)
        specs.append({"version":version,"path":str(p),"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),
                      "source_url":"https://gwosc.org/official/","published_utc":"2026-01-01"})
    tables=load_releases({"schema":"GWTC_CUMULATIVE_MANIFEST_V1","releases":specs})
    rows,summary=reconcile(tables,[{"version_a":"GWTC-4.0","version_b":"GWTC-4.1",
      "from_name":"b","to_name":"renamed","official_identity_evidence_url":"https://gwosc.org/official/"}])
    assert any(r["status"]=="OFFICIAL_RENAME_DECLARED" for r in rows)
    assert any(r["status"]=="REVISED_FIELDS" and r["name_a"]=="a" for r in rows)
    assert all("REMOVED" not in r["status"] for r in rows)
    assert summary[0]["declared_renames"]==1
    specs[0]["sha256"]="0"*64
    with pytest.raises(ValueError,match="SHA"):load_releases({"schema":"GWTC_CUMULATIVE_MANIFEST_V1","releases":specs})
def test_unguarded_rename_rejected():
    with pytest.raises(ValueError,match="official identity evidence"):
        reconcile({"GWTC-4.0":pd.DataFrame(),"GWTC-4.1":pd.DataFrame(),"GWTC-5.0":pd.DataFrame()},
            [{"version_a":"GWTC-4.0","version_b":"GWTC-4.1","from_name":"a","to_name":"b"}])
