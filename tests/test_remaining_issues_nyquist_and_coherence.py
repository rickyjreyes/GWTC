import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from annotate_legacy_nyquist import annotate_legacy
from generate_gwtc4_population_null_matrix import prepare_scan,score_prepared_one_details

def test_archived_legacy_peak_aliased_but_verdict_not_rewritten():
    df=pd.DataFrame({"bin_count":[20,30],"k_best":[32.69747899159664,4.483193277310924],
        "n_star":[23.129202316812606,3.1712745916267093],
        "global_p":[.003,.003],"verdict_label":["PASS","PARTIAL"]})
    a=annotate_legacy(df)
    assert a.nyquist_status.tolist()==["ALIASED_AT_OR_ABOVE_NYQUIST","RESOLVED_BELOW_NYQUIST"]
    assert a.global_p.tolist()==df.global_p.tolist()
    assert a.verdict_label.tolist()==df.verdict_label.tolist()

def test_missing_legacy_inputs_not_promoted():
    a=annotate_legacy(pd.DataFrame({"bin_count":[20],"k_best":[np.nan],"n_star":[np.nan]}))
    assert a.nyquist_status[0]=="UNDETERMINED_MISSING_LEGACY_FIELDS"
    with pytest.raises(ValueError,match="Missing"):
        annotate_legacy(pd.DataFrame({"bin_count":[20]}))

def test_recalculated_historical_observed_details_include_identity_and_n():
    x=np.linspace(1,3,25);y=.2+.1*np.cos(9.7*x)
    df=pd.DataFrame({"event":[str(i) for i in range(25)],
        "log_final_mass_source":x,"final_spin_median":y})
    prep=prepare_scan(df,"log_final_mass_source","final_spin",12,np.linspace(.5,16,100),2)
    q=score_prepared_one_details(prep,y)
    assert len(q["events"].split(","))==25
    assert q["n_obs"]==pytest.approx(q["k_best"]*q["delta_x"]/(2*np.pi))
