import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from gwtc_hierarchical_fixed_mode import hierarchical_statistic,expectation_ratio
def test_constant_mode_cancels_rate_marginalized():
    ell=np.linspace(1,3,100)
    e=[(ell,np.ones(len(ell)))]*2
    r=hierarchical_statistic(e,(ell,np.ones(len(ell))),k=9.7,a=0,b=0,n_generated=100,
         ess_min_event=1,ess_min_injection=1)
    assert r["statistic_2loglr"]==pytest.approx(0,abs=1e-10)
    assert r["numerical_gate"]=="PASS"
def test_selection_correction_can_reduce_naive_evidence():
    e=[(np.full(100,2.0),np.ones(100))]
    inj=(np.full(100,2.0),np.ones(100))
    r=hierarchical_statistic(e,inj,9.7,.175,.343,100,ess_min_event=1,ess_min_injection=1)
    assert r["statistic_2loglr"]==pytest.approx(0,abs=1e-9)
def test_missing_or_negative_weights_fail():
    with pytest.raises(ValueError,match="invalid"):
        expectation_ratio(np.array([1.,2.]),np.array([1.,-1.]),9.7,.1,.2)
def test_low_ess_is_incomplete_not_positive_claim():
    ell=np.array([1.,2.,3.]);w=np.array([100.,0.,0.])
    r=hierarchical_statistic([(ell,w)],(ell,w),9.7,.1,.2,100)
    assert r["numerical_gate"]=="INCOMPLETE"
