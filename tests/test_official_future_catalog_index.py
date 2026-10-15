"""Official-index gate must never auto-authorize a future signal score."""
import json
import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_official_future_catalog_index import audit_index, main


def item(name):
    return dict(name=name, detail_url=f"https://gwosc.org/api/v2/catalogs/{name}?format=json")


def index(*names):
    return dict(next=None, page_number=1, num_pages=1, results_count=len(names),
                results=[item(name) for name in names])


def test_only_old_catalog_is_blocked():
    x=audit_index(index("GWTC-4.1","GWTC-5.0","GWTC"), "hash")
    assert x["first_post_5_0_candidate"] is None
    assert x["status"]=="EXTERNALLY_BLOCKED_NO_POST_5_0_IN_INDEX_SNAPSHOT"
    assert x["signal_scored"] is False and x["event_data_read"] is False


def test_candidate_is_not_authorization():
    x=audit_index(index("GWTC-5.0","GWTC-6.1","GWTC-6.0"), "hash")
    assert x["first_post_5_0_candidate"] == "GWTC-6.0"
    assert x["status"] == "CANDIDATE_REQUIRES_EXTERNAL_RELEASE_ATTESTATION"
    assert x["manifest_prepared"] is False and x["signal_scored"] is False


@pytest.mark.parametrize("bad", [
    dict(next="https://gwosc.org/more", page_number=1,num_pages=2,
         results_count=1,results=[item("GWTC-5.0")]),
    index("GWTC-4.1"),
    index("GWTC-5.0","GWTC-5.0"),
])
def test_incomplete_missing_or_duplicate_index_fails(bad):
    with pytest.raises(ValueError):
        audit_index(bad,"hash")


def test_cli_keeps_snapshot_bytes_and_refuses_overwrite(tmp_path):
    s=tmp_path/"snapshot.json"; before=json.dumps(index("GWTC-5.0")).encode()
    s.write_bytes(before); output=tmp_path/"audit.json"
    assert main(["--snapshot",str(s),"--out",str(output)])==0
    assert s.read_bytes()==before
    assert json.loads(output.read_text())["signal_scored"] is False
    with pytest.raises(SystemExit):
        main(["--snapshot",str(s),"--out",str(output)])
