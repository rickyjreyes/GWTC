import sys
from pathlib import Path
import numpy as np
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from gwtc_hierarchical_baseline_grid import conditional_mixture
from gwtc_posterior_uncertainty import weighted_quantile

def model(name, weight, statistic, mode="frozen", events=None, gate="PASS"):
    return {"model_id":name,"baseline_posterior_weight":weight,"statistic_2loglr":statistic,
        "numerical_gate":gate,"fixed_prediction_sha256":mode,
        "event_names":["one","two"] if events is None else events}
def test_mixture_constant_lr_preserves_lr():
    x=conditional_mixture([model("m1",.2,4),model("m2",.8,4)])
    assert x["conditional_model_averaged_2loglr"]==pytest.approx(4.)
    assert x["calibrated_p"] is None
def test_mixture_not_average_of_log_likelihoods():
    a=conditional_mixture([model("m1",.5,0),model("m2",.5,4)])
    assert a["conditional_model_averaged_2loglr"]==pytest.approx(2*np.log((1+np.exp(2))/2))
    assert a["conditional_model_averaged_2loglr"]!=pytest.approx(2.)
def test_incomplete_and_mismatched_models_rejected():
    with pytest.raises(ValueError,match="Numerically"):
        conditional_mixture([model("m",1.,2,gate="INCOMPLETE")])
    with pytest.raises(ValueError,match="Mismatched"):
        conditional_mixture([model("a",.5,0),model("b",.5,1,events=["other"])])
    with pytest.raises(ValueError,match="Unnormalized"):
        conditional_mixture([model("a",.4,0),model("b",.4,1)])
def test_weighted_posterior_median_uses_weights_not_unweighted_median():
    assert weighted_quantile(np.array([10.,20.,30.]),np.array([.05,.9,.05]),.5)==20.
    with pytest.raises(ValueError):weighted_quantile(np.array([1.]),np.array([-1.]),.5)
