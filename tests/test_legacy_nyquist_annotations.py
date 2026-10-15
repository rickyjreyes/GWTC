"""Issue #15: historical annotations are derived, not retroactive rescans."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from annotate_legacy_nyquist import classify, audit, INPUTS, main


def test_known_historical_20_bin_peak_flagged_aliased():
    r = classify(dict(bin_count="20", k_best="32.69747899159664",
                      n_star="23.129202316812606", verdict_label="PASS"),
                 "historical", 2)
    assert r["sampling_classification"] == "ALIASED_OR_NYQUIST"
    assert float(r["nyquist_k"]) < float(r["historical_k_best"])
    assert r["historical_verdict"] == "PASS"  # retained, not endorsed


def test_historical_30_bin_peak_is_resolved():
    r = classify(dict(bin_count="30", k_best="4.483193277310924",
                      n_star="3.1712745916267093"), "historical", 2)
    assert r["sampling_classification"] == "RESOLVED"


def test_nyquist_boundary_is_not_called_resolved():
    # active width 2*pi => n_star=k; Nyquist=B/2
    r = classify(dict(bin_count="20", k_best="10", n_star="10"), "test", 2)
    assert r["sampling_classification"] == "ALIASED_OR_NYQUIST"


def test_missing_metadata_is_explicit_not_resolved():
    assert classify(dict(bin_count="20", k_best="", n_star=""), "test", 2)[
        "sampling_classification"] == "INSUFFICIENT_METADATA"


def test_archived_inputs_are_untouched(tmp_path):
    root = Path(__file__).resolve().parents[1]
    before = {p: (root / p).read_bytes() for p in INPUTS}
    out = tmp_path / "annotations.csv"
    assert main(["--root", str(root), "--out", str(out)]) == 0
    assert out.exists()
    assert {p: (root / p).read_bytes() for p in INPUTS} == before
    assert any(r["sampling_classification"] == "ALIASED_OR_NYQUIST" for r in audit(root))
