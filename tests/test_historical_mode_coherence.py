import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from historical_mode_coherence import coherence_statistic
from generate_gwtc4_population_null_matrix import prepare_scan,score_prepared_one,score_prepared_one_details,score_prepared_many_details
from gwtc4_wct_subset_scan_compound import harmonic_scan
CFG={"min_jaccard":.5,"k_tolerance":.5,"phase_tolerance_radians":np.pi/4}
def make_row(name,events,k,phase=0):
    return {"scan_id":name,"x_var":"log_final_mass_source","y_var":"final_spin",
     "events":",".join(events),"k_best":k,"phase":phase,"delta_x":2.,
     "n_obs":k*2/(2*np.pi)}
def test_duplicate_membership_not_pseudoreplication():
    a=[str(x) for x in range(12)]
    b=[str(x) for x in range(2,14)]
    df=pd.DataFrame([make_row("a",a,9.7),make_row("a_copy",a,9.7),make_row("b",b,9.8)])
    r=coherence_statistic(df,CFG)
    assert r["deduplicated_scans"]==2 and r["pair_count"]==1 and r["score"]==pytest.approx(1.)
def test_antiphase_is_not_coherent():
    a=[str(x) for x in range(12)];b=[str(x) for x in range(2,14)]
    df=pd.DataFrame([make_row("a",a,9.7),make_row("b",b,9.8,np.pi)])
    assert coherence_statistic(df,CFG)["score"]==0
def test_vectorized_peak_details_match_original_harmonic_fit():
    x=np.linspace(1,3,30); y=.5+.2*np.cos(9.7*x-.3)+.01*x
    df=pd.DataFrame({"event":[str(i) for i in range(30)],
                     "log_final_mass_source":x,"final_spin_median":y})
    kg=np.linspace(.5,16,100)
    prep=prepare_scan(df,"log_final_mass_source","final_spin",12,kg,2)
    oldk,oldd=score_prepared_one(prep,y)
    details=score_prepared_one_details(prep,y)
    many=score_prepared_many_details(prep,np.stack([y,y]))
    _,historic,_,_=harmonic_scan(x,y,None,kg,2)
    assert details["k_best"]==pytest.approx(oldk)
    assert details["delta_chi2"]==pytest.approx(oldd,abs=1e-6)
    assert details["phase"]==pytest.approx(historic["phase"],abs=1e-5)
    assert np.allclose(many[:,0],oldk) and np.allclose(many[:,1],oldd,atol=1e-6)
