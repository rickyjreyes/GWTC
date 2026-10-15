import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from evaluate_future_gwtc import release_ledger

def test_complete_new_identity_dispositions_before_scoring():
    future=pd.DataFrame({"commonName":["old","valid","low","badmass","nonprimary"],
      "M_chirp":[12,15,18,-1,99],"p_astro":[.9,.6,.1,.8,.9],
      "is_primary_entry":[True,True,True,True,False]})
    h=pd.DataFrame({"commonName":["old","absent_old"]})
    d=release_ledger(future,h,{"selection":{"p_astro_min_inclusive":.5}})
    assert d["eligible_names"]==["badmass","valid"]
    assert d["eligible_invalid_mass_names"]==["badmass"]
    assert d["new_names_below_or_missing_pastro"]==["low"]
    assert d["historical_names_not_listed_as_primary"]==["absent_old"]
    assert d["n_primary_rows"]==4
def test_duplicate_primary_fails_closed():
    f=pd.DataFrame({"commonName":["x","x"],"M_chirp":[10,12],"p_astro":[.9,.8],"is_primary_entry":[True,True]})
    with pytest.raises(ValueError,match="Duplicate"):release_ledger(f,pd.DataFrame({"commonName":["old"]}),{"selection":{"p_astro_min_inclusive":.5}})
