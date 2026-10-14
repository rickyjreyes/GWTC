import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from evaluate_future_gwtc import selected_names

CFG={"selection":{"p_astro_min_inclusive":0.5}}
HIST=pd.DataFrame({"commonName":["old","old2"]})
def test_only_new_primary_and_frozen_cut():
    future=pd.DataFrame({"commonName":["old","new","bad","new2","low"],
       "M_chirp":[9,10,11,12,13],"p_astro":[.99,.5,.9,.9,.49],
       "is_primary_entry":[True,True,False,True,True]})
    assert selected_names(future,HIST,CFG).commonName.tolist()==["new","new2"]
def test_duplicate_primary_fails_before_scoring():
    f=pd.DataFrame({"commonName":["x","x"],"M_chirp":[10,12],"p_astro":[.9,.9],
                    "is_primary_entry":[True,True]})
    with pytest.raises(ValueError,match="Duplicate"):selected_names(f,HIST,CFG)
