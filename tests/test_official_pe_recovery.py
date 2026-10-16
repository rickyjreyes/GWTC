import copy,hashlib,importlib.util,json,sys
from pathlib import Path
import h5py,numpy as np,pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from recover_official_pe_inputs import preferred_event,verify_file
from run_official_pe_shard import project,accepted_projection

def official_fixture():
 return {'events':{'GW240420_175625-v1':{'commonName':'GW240420_175625','catalog.shortName':'GWTC-5.0','parameters':{
  'preferred':{'pipeline_type':'pe','is_preferred':True,'waveform_family':'C00:IMRPhenomXPHM-SpinTaylor',
   'data_url':'https://zenodo.org/api/records/20348005/files/IGWN-GWTC5p0-29ebe06b7_25-GW240420_175625-combined_PEDataRelease.hdf5/content',
   'chirp_mass_source_unit':'M_sun','chirp_mass_source':24.7}}}}}

def test_units_identity_host_and_preferred_ambiguity_fail_closed():
 d=official_fixture();name='GW240420_175625'
 assert preferred_event(d,name)['frame']=='source'
 for field,value in [('chirp_mass_source_unit','kg'),('data_url','https://example.com/GW240420_175625/file')]:
  bad=copy.deepcopy(d);bad['events'][name+'-v1']['parameters']['preferred'][field]=value
  with pytest.raises(ValueError):preferred_event(bad,name)
 bad=copy.deepcopy(d);bad['events'][name+'-v1']['parameters']['second']=bad['events'][name+'-v1']['parameters']['preferred'].copy()
 with pytest.raises(ValueError):preferred_event(bad,name)

def test_lossless_projection_and_corrupt_raw_checksum(tmp_path):
 raw=tmp_path/'raw.hdf5';values=np.zeros(4,dtype=[('chirp_mass_source','f8'),('peak_luminosity','f8')]);values['chirp_mass_source']=[1.,2.,3.,4.]
 with h5py.File(raw,'w') as f:f.create_dataset('C00:IMRPhenomXPHM-SpinTaylor/posterior_samples',data=values)
 row=preferred_event(official_fixture(),'GW240420_175625');row.update(bytes=raw.stat().st_size,official_md5=hashlib.md5(raw.read_bytes()).hexdigest())
 receipt=project(row,raw,tmp_path)
 assert receipt['posterior_n']==receipt['weight_ess']==4
 assert receipt['projection_lossless_verified'] and receipt['secondary_analysis_ready'] is False
 assert receipt['secondary_fields']['peak_luminosity']['status']=='PRESENT_EXACT_NAME'
 assert accepted_projection(row,receipt,tmp_path)==receipt
 with pytest.raises(ValueError,match='identity mismatch'):accepted_projection({**row,'waveform':'wrong'},receipt,tmp_path)
 with pytest.raises(ValueError,match='checksum mismatch'):accepted_projection(row,{**receipt,'projection_sha256':'0'*64},tmp_path)
 with h5py.File(tmp_path/'GW240420_175625.hdf5') as f:np.testing.assert_array_equal(f['posterior_samples']['chirp_mass_source'],values['chirp_mass_source'])
 bad={**row,'official_md5':'0'*32}
 with pytest.raises(ValueError,match='MD5'):verify_file(bad,raw)

def test_detector_frame_only_samples_are_rejected(tmp_path):
 raw=tmp_path/'wrong.hdf5'
 with h5py.File(raw,'w') as f:f.create_dataset('C00:IMRPhenomXPHM-SpinTaylor/posterior_samples',data=np.zeros(3,dtype=[('chirp_mass','f8')]))
 row=preferred_event(official_fixture(),'GW240420_175625');row.update(bytes=raw.stat().st_size,official_md5=hashlib.md5(raw.read_bytes()).hexdigest())
 with pytest.raises(ValueError,match='source-frame'):verify_file(row,raw)
