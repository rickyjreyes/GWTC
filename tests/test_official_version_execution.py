import copy,sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from gwtc_official_version_execution import normalize,reconcile

def fixture(name='GW240101_010101',version='GWTC-4.0'):
 return {'events':{name+'-v1':{'commonName':name,'version':1,'catalog.shortName':version,'GPS':123.,
  'jsonurl':'https://gwosc.org/eventapi/json/'+version+'/'+name+'/v1','chirp_mass_source':20.,'chirp_mass_source_unit':'M_sun','p_astro':.99}}}

def test_invalid_units_and_duplicate_physical_identity_rejected():
 obj=fixture();next(iter(obj['events'].values()))['chirp_mass_source_unit']='kg'
 with pytest.raises(ValueError,match='units'):normalize(obj,'GWTC-4.0')
 obj=fixture();obj['events'].update(fixture('GW240101_010102')['events'])
 with pytest.raises(ValueError,match='same-GPS'):normalize(obj,'GWTC-4.0')

def test_missing_mass_retained_and_rename_not_silently_merged():
 a=normalize(fixture(),'GWTC-4.0');obj=fixture('GW240101_010102','GWTC-4.1');next(iter(obj['events'].values()))['chirp_mass_source']=None
 b=normalize(obj,'GWTC-4.1');r=reconcile(a,b)
 assert len(b)==1 and b.M_chirp.isna().all()
 assert len(r['same_gps_rename_candidates'])==1 and r['automatic_alias_merges']==0 and not r['physical_removal_inferred']

def test_shared_name_changed_physical_identity_rejected():
 a=normalize(fixture(),'GWTC-4.0');b=normalize(fixture(version='GWTC-4.1'),'GWTC-4.1');b.loc[0,'GPS']=124.
 with pytest.raises(ValueError,match='GPS identity'):reconcile(a,b)
