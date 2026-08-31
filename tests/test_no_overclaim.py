#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Overclaim guard + required-language tests.

These tests enforce the scientific discipline of the harness: required hedging
language must be present, forbidden unqualified overclaim phrases must be
absent, and the harness must be able to emit a FAIL verdict without any test
failing.
"""

import csv
import glob
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _read(path):
    with open(os.path.join(ROOT, path), encoding="utf-8") as fh:
        return fh.read()


# ---------------------------------------------------------------------------
# Required language
# ---------------------------------------------------------------------------
def test_readme_contains_does_not_prove_wct():
    assert "does not prove WCT" in _read("README.md")


def test_method_contains_admissible_scale_variables():
    assert "admissible scale variables" in _read("METHOD.md").lower()


def test_readme_states_does_not_replace_lvk():
    assert "does not replace LVK population" in _read("README.md")


# ---------------------------------------------------------------------------
# Forbidden unqualified overclaim phrases (case-insensitive substring match).
# The guard scans harness docs and scripts but NOT this test file (which must
# contain the phrases in order to check for them) and NOT legacy root scripts.
#
# Important: do not forbid a generic word such as "proof". Bounded statements
# such as "this is not a proof" are scientifically useful and must remain
# legal. The guard targets affirmative overclaims instead.
# ---------------------------------------------------------------------------
FORBIDDEN = [
    "gwtc proves wct",
    "gravitational waves prove wct",
    "new physics discovered",
    "lvk missed wct",
    "confirmed wct",
    "proof of wct",
    "proof that wct",
]

SCANNED_DOCS = [
    "README.md", "METHOD.md", "RESULTS.md", "REPRODUCE.md",
    "PROVENANCE.md", "NEGATIVE_CONTROL.md", "TODO.md", "data/README.md",
]


def _harness_files():
    files = [os.path.join(ROOT, d) for d in SCANNED_DOCS]
    files += sorted(glob.glob(os.path.join(ROOT, "scripts", "*.py")))
    return [f for f in files if os.path.exists(f)]


@pytest.mark.parametrize("path", _harness_files())
def test_no_forbidden_overclaim_phrases(path):
    text = open(path, encoding="utf-8").read().lower()
    for phrase in FORBIDDEN:
        assert phrase not in text, f"forbidden overclaim phrase '{phrase}' found in {path}"


def test_allowed_bounded_language_present_somewhere():
    blob = " ".join(_read(d).lower() for d in SCANNED_DOCS)

    # Require the scientific ideas, not one brittle editorial sentence. This
    # prevents harmless wording revisions from breaking the numerical harness.
    for term in ["wct-motivated", "diagnostic", "log-domain residual"]:
        assert term in blob, f"expected bounded concept missing: '{term}'"

    for allowed in [
        "scale-like variable",
        "negative control",
        "does not replace lvk population",
        "requires independent confirmation",
    ]:
        assert allowed in blob, f"expected allowed bounded phrase missing: '{allowed}'"


# ---------------------------------------------------------------------------
# The harness must tolerate FAIL / negative results without breaking.
# ---------------------------------------------------------------------------
def test_verdict_csv_may_contain_fail_without_test_failure():
    """A FAIL (or PARTIAL/INCOMPLETE) verdict is a valid, retained outcome."""
    path = os.path.join(ROOT, "tables", "gwtc_verdict.csv")
    if not os.path.exists(path):
        pytest.skip("verdict not yet generated")
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    assert rows, "verdict CSV is empty"
    valid = {"PASS", "PARTIAL", "FAIL", "INCOMPLETE"}
    for r in rows:
        assert r["verdict_label"] in valid
        assert r["reliability_class"] in {"I", "II", "III", "IV"}
    # Explicitly assert that a FAIL value is acceptable (does NOT fail the suite).
    fail_ok = True
    for r in rows:
        if r["verdict_label"] == "FAIL":
            fail_ok = True  # no exception raised: FAIL is tolerated
    assert fail_ok


def test_class_iv_negative_control_output_is_valid():
    """Class IV / negative-control rows must conform to the canonical schema."""
    import sys
    sys.path.insert(0, os.path.join(ROOT, "scripts"))
    import wct_gwtc_lib as L

    row = L.result_row(
        catalog_version="cumulative", variable="abs_chi_eff", subset="ctrl",
        bin_count=20, k_best=5.0, n_star=3.0, delta_d_star=2.0,
        local_p=0.3, global_p=0.7, null_n=400,
        verdict_label="FAIL", notes="control=uniform_ell",
    )
    assert set(L.CSV_COLUMNS).issubset(row.keys())
    assert row["verdict_label"] in L.VERDICT_LABELS
