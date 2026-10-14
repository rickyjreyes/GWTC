import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from gwtc_catalog_version_audit import version_index,version_diff,FIELDS
def row(ver,name,m):
    x={"catalog_version":ver,"commonName":name}
    for f in FIELDS:x[f]=float("nan")
    x["M_chirp"]=m
    return x
def test_version_diff_never_calls_missing_event_removed():
    df=pd.DataFrame([row("GWTC-4.0","a",10),row("GWTC-4.0","b",20),
                     row("GWTC-4.1","a",11),row("GWTC-4.1","c",30)])
    a=version_index(df,"GWTC-4.0");b=version_index(df,"GWTC-4.1")
    d={r["event_name"]:r for r in version_diff(a,b,"GWTC-4.0","GWTC-4.1")}
    assert d["a"]["status"]=="revised" and d["a"]["changed_fields"]=="M_chirp"
    assert d["b"]["status"]=="not_listed_in_b"
    assert d["c"]["status"]=="not_listed_in_a"
def test_duplicate_version_identity_fails():
    df=pd.DataFrame([row("GWTC-4.0","a",10),row("GWTC-4.0","a",11)])
    with pytest.raises(ValueError,match="Ambiguous"):
        version_index(df,"GWTC-4.0")
