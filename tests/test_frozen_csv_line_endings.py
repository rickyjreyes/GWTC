import hashlib
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from gwtc_source_provenance import hashes,validate_frozen_csv

def test_frozen_newline_equivalence_is_explicit_not_silent(tmp_path):
    p=tmp_path/"original.csv"
    crlf=b"identity,mass\r\na,10\r\nb,20\r\n"
    p.write_bytes(crlf)
    expected=hashlib.sha256(crlf).hexdigest()
    exact=validate_frozen_csv(p,expected)
    assert exact["verification_status"]=="EXACT_ORIGINAL_BYTES"
    p.write_bytes(crlf.replace(b"\r\n",b"\n"))
    eol=validate_frozen_csv(p,expected)
    assert eol["verification_status"]=="DECLARED_HASH_MATCHES_ONLY_LF_CRLF_TRANSFORM"
    assert eol["checkout_sha256"]!=expected
    assert eol["normalized_crlf_sha256"]==expected
    p.write_bytes(p.read_bytes().replace(b"b,20",b"b,21"))
    with pytest.raises(ValueError,match="beyond permitted newline"):
        validate_frozen_csv(p,expected)
def test_rejects_bare_cr_format(tmp_path):
    p=tmp_path/"bad.csv";p.write_bytes(b"a\rb")
    with pytest.raises(ValueError,match="bare-CR"):hashes(p)
