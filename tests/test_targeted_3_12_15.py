import copy
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from make_publication_figures import validated_v2_nyquist
from plot_legacy_nyquist import validate_overlay,make_figure
from gwtc_baseline_v2_acceptance import archived_consistency,baseline_gate
from gwtc_version_score_decomposition import scored_versions,decompose

def test_real_frozen_v2_presentation_is_resolved():
    v2=json.loads((ROOT/"tables/gwtc_frozen_mode.json").read_text())
    assert v2["bins"]==40
    assert v2["k_star"]<validated_v2_nyquist(v2)
    changed=copy.deepcopy(v2)
    changed["k_star"]=changed["nyquist_k"]*1.1
    with pytest.raises(ValueError,match="Aliased"):
        validated_v2_nyquist(changed)

def test_readonly_legacy_figure_refuses_false_alias_or_overwrite(tmp_path):
    frame=pd.read_csv(ROOT/"results/legacy_nyquist/gwtc_bin_stress_annotated.csv")
    checked=validate_overlay(frame)
    assert checked.loc[checked.bin_count==20,"nyquist_status"].iloc[0].startswith("ALIASED")
    bad=frame.copy();bad.loc[bad.bin_count==20,"nyquist_status"]="RESOLVED_BELOW_NYQUIST"
    with pytest.raises(ValueError,match="Contradictory"):
        validate_overlay(bad)
    target=tmp_path/"historical_audit.svg"
    make_figure(frame,target)
    assert "ALIASED" in target.read_text()
    with pytest.raises(FileExistsError):make_figure(frame,target)

def test_real_archived_baseline_v2_provenance_checks():
    audit=archived_consistency(ROOT)
    assert audit["status"]=="ARCHIVED_CONSISTENT_NOT_PROSPECTIVE_REPLICATION"
    assert audit["v2_selected_polynomial_degree"]==7
    assert audit["holdout_positive_n"]==104
    assert audit["v2_frozen_k"]<audit["v2_nyquist_k"]
    gate=baseline_gate(root=ROOT)
    assert gate["classification"]=="PARTIAL_CLASS_III_NOT_PROMOTED"
    assert gate["gates"]["posterior_measurement_uncertainty"]["status"]=="NOT_RUN_NO_OFFICIAL_POSTERIOR_DIAGNOSTIC"

def test_tagged_score_difference_is_additive_even_when_selection_flips():
    rows=[]
    for ver,records in {
       "GWTC-4.0":[("a",10,.99),("b",20,.99)],
       "GWTC-4.1":[("a",12,.99),("b",20,.1),("c",30,.99)],
       "GWTC-5.0":[("a",12,.99),("c",30,.99),("d",40,.99)]}.items():
        for name,m,p in records:rows.append({"catalog_version":ver,"commonName":name,"M_chirp":m,"p_astro":p})
    mode={"k":9.7,"a_cos":.175,"b_sin":.343,"normalizer_Z":1.05}
    versions=scored_versions(pd.DataFrame(rows),mode)
    ledger,pairs=decompose(versions)
    for x in pairs:
        assert x["score_change_b_minus_a"]==pytest.approx(sum(x["decomposition"].values()),abs=1e-10)
    assert pairs[0]["shared_selection_flips"]==1
    assert any(x["identity_category"]=="not_listed_in_a" for x in ledger)
    assert any(x["identity_category"]=="not_listed_in_b" for x in ledger)
    assert pairs[0]["score_b"]==pytest.approx(float(versions["GWTC-4.1"]["score"].sum()))
