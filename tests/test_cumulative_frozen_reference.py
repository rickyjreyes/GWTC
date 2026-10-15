import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
import gwtc_cumulative_frozen_reference as M

def test_frozen_cumulative_reference_keeps_versions_separate(monkeypatch):
    monkeypatch.setattr(M,"fixed_model_null",lambda **kw: np.full(kw["null_n"],-1000.))
    tables={
      "GWTC-4.0":pd.DataFrame({"M_chirp":[10.,-1.,30.],"p_astro":[.8,.8,.1]}),
      "GWTC-4.1":pd.DataFrame({"M_chirp":[11.,12.,13.],"p_astro":[.8,.8,.8]}),
      "GWTC-5.0":pd.DataFrame({"M_chirp":[0.,22.],"p_astro":[.8,.8]})}
    v3={"train_ell":[2.,3.],"selected_bandwidth":.5}
    v4={"fixed_mode":{"k":9.7,"a_cos":.175,"b_sin":.343,"normalizer_Z":1.05}}
    rows=M.evaluate_cumulative(tables,v3,v4,1000,42)
    assert [r["positive_pastro_and_mass_n"] for r in rows]==[1,3,1]
    assert [r["missing_or_nonpositive_mass_after_pastro_cut"] for r in rows]==[1,0,1]
    assert all(r["frozen_k"]==9.7 for r in rows)
    assert all(r["fixed_reference_p"]==pytest.approx(1/1001) for r in rows)
    assert [r["seed"] for r in rows]==[42,1000045,2000048]
    assert all(r["population_or_selection_corrected"] is False and r["prospective"] is False for r in rows)
    with pytest.raises(ValueError,match="1000"):
        M.evaluate_cumulative(tables,v3,v4,999,42)

def test_missing_valid_mass_reports_no_reference_p(monkeypatch):
    def forbidden(**kw):raise AssertionError("No null on empty sample")
    monkeypatch.setattr(M,"fixed_model_null",forbidden)
    tables={v:pd.DataFrame({"M_chirp":[np.nan,-1.],"p_astro":[.9,.99]}) for v in M.VERSIONS}
    rows=M.evaluate_cumulative(tables,{"train_ell":[1.],"selected_bandwidth":.5},
      {"fixed_mode":{"k":9.7,"a_cos":.1,"b_sin":.2,"normalizer_Z":1.1}},1000,3)
    assert all(r["fixed_reference_p"] is None and r["statistic_2loglr"] is None for r in rows)
