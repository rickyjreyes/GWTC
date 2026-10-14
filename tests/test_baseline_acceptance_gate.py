import sys
import json
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from gwtc_baseline_v2_acceptance import baseline_gate,sha
def test_no_data_preserves_partial_and_missing_gate(tmp_path):
    x=baseline_gate(root=tmp_path)
    assert x["classification"]=="PARTIAL_CLASS_III_NOT_PROMOTED"
    assert not x["promotion"]
    assert x["gates"]["posterior_measurement_uncertainty"]["status"]=="NOT_RUN_NO_OFFICIAL_POSTERIOR_DIAGNOSTIC"
def test_posterior_sha_and_frozen_binding_gate(tmp_path):
    folder=tmp_path/"tables";folder.mkdir()
    frozen=folder/"gwtc_v4_frozen_structured_null.json";frozen.write_text("{}")
    path=tmp_path/"out.json";path.write_text(json.dumps({"schema":"GWTC_V4_POSTERIOR_SENSITIVITY_V1",
      "frozen_sha256":sha(frozen),"coverage":"COMPLETE","input_events":2,
      "holdout_expected_positive_mass":2,"missing_event_names":[],"calibrated_p":None}))
    q=baseline_gate(path,sha(path),root=tmp_path)
    assert q["scientific_verdict"]=="POSTERIOR_DIAGNOSTIC_DONE_REQUIRE_SCIENTIFIC_INTERPRETATION_AND_ORIGINAL_GATE_AUDIT"
    assert q["classification"]=="PARTIAL_CLASS_III_NOT_PROMOTED"
    with pytest.raises(ValueError,match="SHA"):
        baseline_gate(path,"0"*64,root=tmp_path)
