import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from gwtc_secondary_observable_preflight import preflight
def test_exact_field_without_documented_units_stays_blocked():
    data=pd.DataFrame([{"event_name":"event","observable":"E_rad","status":"PRESENT_EXACT_NAME",
       "field_names":"radiated_energy","file_sha256":"fakehash"}])
    defs={"schema":"GWTC_SECONDARY_OFFICIAL_DEFINITIONS_V1","observables":[{"observable":"E_rad","exact_field":"radiated_energy"}]}
    r=preflight(data,defs)
    assert r.preflight_status.iloc[0]=="BLOCKED"
    assert "MISSING_UNITS" in r.reasons.iloc[0]
def test_only_explicit_reviewed_provenance_is_conditional_not_final():
    data=pd.DataFrame([{"event_name":"event","observable":"E_rad","status":"PRESENT_EXACT_NAME",
       "field_names":"radiated_energy","file_sha256":"fakehash"}])
    defs={"schema":"GWTC_SECONDARY_OFFICIAL_DEFINITIONS_V1","observables":[{"observable":"E_rad",
       "exact_field":"radiated_energy","official_definition_url":"https://gwosc.org/x",
       "release_id":"Official Release","definition_sha256":"a"*64,"units":"declared",
       "frame_convention":"declared","independent_definition_reviewed":True}]}
    r=preflight(data,defs)
    assert r.preflight_status.iloc[0]=="CONDITIONAL_DEFINITION_REVIEWED_NOT_SCIENTIFIC_RESULT"
    data.loc[0,"status"]="ABSENT"
    assert preflight(data,defs).preflight_status.iloc[0]=="BLOCKED"
