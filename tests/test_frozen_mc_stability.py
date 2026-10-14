import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from run_frozen_mc_stability import clopper_pearson
def test_binomial_edge_intervals():
    low,up=clopper_pearson(0,100)
    assert low==0 and 0<up<.1
    low,up=clopper_pearson(100,100)
    assert up==1 and .9<low<1
def test_nonzero_binomial_interval():
    lo,hi=clopper_pearson(55,10000)
    assert lo<.0055<hi
