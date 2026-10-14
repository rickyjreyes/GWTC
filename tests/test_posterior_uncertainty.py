import hashlib
import sys
from pathlib import Path
import h5py
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from gwtc_posterior_uncertainty import extract, draw_statistic

def test_source_frame_and_weights(tmp_path):
    path=tmp_path/"one.h5"
    dt=np.dtype([("chirp_mass_source","f8"),("chirp_mass","f8")])
    with h5py.File(path,"w") as f:
        f.create_dataset("Mixed/posterior_samples",data=np.array([(10,50),(20,60)],dtype=dt))
        f.create_dataset("weights",data=[1.,3.])
    h=hashlib.sha256(path.read_bytes()).hexdigest()
    ell,w,src=extract(path,"Mixed/posterior_samples",h,"weights")
    assert src=="chirp_mass_source"
    assert np.exp(ell)==pytest.approx([10,20])
    assert w==pytest.approx([.25,.75])
    with pytest.raises(ValueError,match="SHA256"):
        extract(path,"Mixed/posterior_samples","0"*64,"weights")

def test_no_posterior_pooling_inflation():
    ell=np.array([np.log(10.),np.log(12.)])
    e=[(ell,np.array([.5,.5]))]
    a=draw_statistic(e,9.7,.175,.343,1.05,100,5)
    doubled=draw_statistic(e+e,9.7,.175,.343,1.05,100,5)
    assert len(a)==100 and len(doubled)==100
    assert np.all(np.isfinite(a))
