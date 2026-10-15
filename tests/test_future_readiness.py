"""Issue #18 preflight: no science scoring, preserve external release gate."""
import sys
from pathlib import Path
import pandas as pd
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from check_future_gwtc_readiness import check_candidate, main


def test_blocked_status_is_not_replication(capsys):
    assert main(["--config", str(Path(__file__).resolve().parents[1] /
                                 "protocols/future_gwtc_v1.json")]) == 0
    printed = capsys.readouterr().out
    assert "EXTERNALLY_BLOCKED_NO_ATTESTED_UNSEEN_RELEASE" in printed
    assert '"scores_computed": false' in printed


def test_old_release_rejected_before_access(tmp_path):
    with pytest.raises(ValueError, match="historical"):
        check_candidate(tmp_path / "not-needed.csv", {}, "GWTC-5.0",
                        "https://gwosc.org/GWTC-5.0/", "2026-01-01")


def test_no_official_provenance_fails(tmp_path):
    with pytest.raises(ValueError, match="Official source"):
        check_candidate(tmp_path / "not-needed.csv", {}, "GWTC-X",
                        "", "2026-12-01")


def test_candidate_only_counts_new_primary_without_evaluation(tmp_path, monkeypatch):
    import check_future_gwtc_readiness as mod
    monkeypatch.setattr(mod, "check_historical", lambda cfg: None)
    monkeypatch.setattr(mod, "sha", lambda path: "fake-checksum")
    cfg = dict(schema="GWTC_FUTURE_FROZEN_PROTOCOL_V1",
               historical_table_path="history.csv",
               selection=dict(p_astro_min_inclusive=.5))
    monkeypatch.setattr(mod.pd, "read_csv",
                        lambda p: pd.DataFrame(dict(commonName=["old"])) if
                        str(p) == "history.csv" else pd.DataFrame(dict(
                            commonName=["old", "new", "low"],
                            M_chirp=[10., 12., 14.],
                            p_astro=[.99, .6, .4],
                            is_primary_entry=[True, True, True])))
    obj = check_candidate(tmp_path / "candidate.csv", cfg, "GWTC-FUTURE",
                          "https://gwosc.org/future", "2027-01-01")
    assert obj["candidate_primary_eligible_names"] == ["new"]
    assert obj["scores_computed"] is False
    assert obj["holdout_evaluated"] is False
    assert obj["status"].startswith("CANDIDATE_PENDING")
