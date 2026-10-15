"""Scheduled official GWTC index watch remains strictly non-scoring."""
import hashlib
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from watch_future_gwtc_release import screen_snapshot, write_snapshot, main


def item(name):
    return {"name": name, "detail_url": f"https://gwosc.org/api/v2/catalogs/{name}?format=json"}


def snapshot(*names):
    return json.dumps(dict(next=None, page_number=1, num_pages=1,
                           results_count=len(names), results=[item(n) for n in names])).encode()


def test_blocked_snapshot_is_archived_byte_exact(tmp_path):
    raw = snapshot("GWTC-4.1", "GWTC-5.0")
    target = tmp_path / "watch"
    report = write_snapshot(target, raw, observed_utc="2026-10-15T20:00:00+00:00")
    assert (target / "official_gwosc_index_raw.json").read_bytes() == raw
    assert report["snapshot_sha256"] == hashlib.sha256(raw).hexdigest()
    assert report["status"] == "EXTERNALLY_BLOCKED_NO_POST_5_0_IN_INDEX_SNAPSHOT"
    assert report["signal_scored"] is False and report["event_data_read"] is False
    assert report["future_result"] is None and report["first_release_attested"] is False
    assert report["observed_utc_not_publication_time"] == "2026-10-15T20:00:00+00:00"
    with pytest.raises(FileExistsError):
        write_snapshot(target, raw, observed_utc="later")


def test_candidate_is_only_alert_not_prepare(tmp_path):
    raw = snapshot("GWTC-5.0", "GWTC-6.0")
    r = screen_snapshot(raw, observed_utc="2026-10-15T20:00:00+00:00")
    assert r["first_post_5_0_candidate"] == "GWTC-6.0"
    assert r["workflow_action"] == "EXTERNAL_CHRONOLOGY_AND_EVENT_IDENTITY_REVIEW_REQUIRED"
    assert r["manifest_prepared"] is False and r["signal_scored"] is False
    assert r["future_manifest_sha256"] is None


def test_malformed_and_pagination_fail_without_archiving(tmp_path):
    for raw in (b"{}", b"not json", json.dumps(dict(next="more",page_number=1,
                   num_pages=2,results_count=1,results=[item("GWTC-5.0")])).encode()):
        target = tmp_path / "incomplete"
        with pytest.raises((ValueError, json.JSONDecodeError)):
            write_snapshot(target, raw, observed_utc="now")
        assert not target.exists()


def test_cli_offline_snapshot_no_overwrite(tmp_path):
    raw = snapshot("GWTC-5.0")
    source = tmp_path / "source.json"; source.write_bytes(raw)
    out = tmp_path / "fresh"
    assert main(["--snapshot", str(source), "--out-dir", str(out)]) == 0
    assert source.read_bytes() == raw
    with pytest.raises(SystemExit):
        main(["--snapshot", str(source), "--out-dir", str(out)])
