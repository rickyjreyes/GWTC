import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from audit_official_gwtc_observables import audit_field_names
def test_fields_not_inferred_from_mass_or_other_aliases():
    f=audit_field_names(["final_mass_source","mass_1_source","mass_2_source"])
    assert f["E_rad"]["status"]=="ABSENT"
    assert f["peak_luminosity"]["status"]=="ABSENT"
def test_exact_official_field_presence():
    f=audit_field_names(["radiated_energy","peak_luminosity"])
    assert f["E_rad"]["status"]=="PRESENT_EXACT_NAME"
    assert f["peak_luminosity"]["status"]=="PRESENT_EXACT_NAME"
