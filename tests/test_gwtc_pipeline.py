#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pipeline correctness tests for the WCT-motivated GWTC diagnostic."""

import math
import os
import subprocess
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)

import wct_gwtc_lib as L  # noqa: E402


# ---------------------------------------------------------------------------
# p_global formula
# ---------------------------------------------------------------------------
def test_global_p_formula_basic():
    # T_null all below T_obs -> p = 1/(1+N)
    t_null = np.zeros(99)
    assert L.global_p_value(t_null, t_obs=10.0) == pytest.approx(1 / 100)


def test_global_p_formula_counts_ties_as_ge():
    # 4 of 9 nulls >= T_obs -> (1+4)/(1+9) = 0.5
    t_null = np.array([1, 2, 3, 4, 5, 5, 5, 5, 5], dtype=float)
    assert L.global_p_value(t_null, t_obs=5.0) == pytest.approx((1 + 5) / (1 + 9))


def test_global_p_matches_run_poisson_null_definition():
    # The null routine must implement exactly (1 + #{T_null>=T_obs})/(1+N).
    rng = np.random.default_rng(0)
    z = np.exp(rng.normal(2.0, 0.5, size=200))
    field = L.build_density_field(z, bin_count=20)
    k_grid = L.make_k_grid(0.5, 20.0, 30)
    scan = L.scan_k(field, k_grid)
    res = L.run_poisson_null(field, k_grid, scan.delta_d_star, null_n=50, seed=1)
    expected = (1 + int(np.sum(res.t_null >= scan.delta_d_star))) / (1 + 50)
    assert res.global_p == pytest.approx(expected)
    assert 0.0 < res.global_p <= 1.0


# ---------------------------------------------------------------------------
# n_star formula
# ---------------------------------------------------------------------------
def test_n_star_formula():
    assert L.compute_n_star(k_star=2 * math.pi, delta_ell_active=1.0) == pytest.approx(1.0)
    assert L.compute_n_star(k_star=4 * math.pi, delta_ell_active=2.0) == pytest.approx(4.0)


def test_scan_n_star_consistent_with_formula():
    rng = np.random.default_rng(2)
    z = np.exp(rng.uniform(0.0, 3.0, size=300))
    field = L.build_density_field(z, bin_count=25)
    scan = L.scan_k(field, L.make_k_grid(0.5, 30.0, 40))
    assert scan.n_star == pytest.approx(
        scan.k_star * field.delta_ell_active / (2 * math.pi))


# ---------------------------------------------------------------------------
# deviance / GLM sanity
# ---------------------------------------------------------------------------
def test_poisson_deviance_zero_at_truth():
    y = np.array([3.0, 5.0, 2.0, 8.0])
    assert L.poisson_deviance(y, y) == pytest.approx(0.0, abs=1e-9)


def test_delta_d_nonnegative():
    rng = np.random.default_rng(3)
    z = np.exp(rng.normal(1.5, 0.4, size=250))
    field = L.build_density_field(z, bin_count=20)
    mu0 = L.fit_baseline(field)
    for k in [1.0, 5.0, 12.0]:
        dd = L.delta_d_at_k(field.counts, field.centers, mu0, k)
        assert dd >= -1e-6  # adding a mode cannot worsen the in-sample fit


def test_density_field_active_support():
    z = np.array([math.e ** 0.0, math.e ** 1.0, math.e ** 2.0])
    field = L.build_density_field(z, bin_count=4)
    assert field.delta_ell_active == pytest.approx(2.0)
    assert field.ell_min == pytest.approx(0.0)
    assert field.ell_max == pytest.approx(2.0)


def test_build_density_field_rejects_nonpositive():
    z = np.array([-1.0, 0.0, np.nan, np.inf])
    with pytest.raises(ValueError):
        L.build_density_field(z, bin_count=5)


# ---------------------------------------------------------------------------
# analyze_cell verdict labels (incl. valid Class IV / negative-control output)
# ---------------------------------------------------------------------------
def test_analyze_cell_incomplete_on_few_events():
    row = L.analyze_cell([2.0, 3.0, 4.0], variable="M_chirp",
                         catalog_version="x", subset="s", bin_count=10,
                         k_grid=L.make_k_grid(0.5, 10, 5), null_n=10, min_events=20)
    assert row["verdict_label"] == "INCOMPLETE"
    # canonical schema present
    for col in L.CSV_COLUMNS:
        assert col in row


def test_per_result_verdict_allows_fail_and_incomplete():
    assert L.per_result_verdict(0.9, null_n=1000) == "FAIL"
    assert L.per_result_verdict(None, null_n=0, have_data=False) == "INCOMPLETE"
    assert L.per_result_verdict(0.001, null_n=1000) == "PASS"
    assert L.per_result_verdict(0.001, null_n=50) == "PARTIAL"  # too few nulls


def test_bounded_variable_flagged_diagnostic_only():
    rng = np.random.default_rng(5)
    vals = np.abs(rng.normal(0, 0.3, size=200))
    vals = vals[vals > 0]
    row = L.analyze_cell(vals, variable="abs_chi_eff", catalog_version="x",
                         subset="s", bin_count=20, k_grid=L.make_k_grid(0.5, 20, 20),
                         null_n=10)
    assert "bounded_non_scale_diagnostic_only" in row["notes"]


# ---------------------------------------------------------------------------
# scripts expose a working --help
# ---------------------------------------------------------------------------
SCRIPT_NAMES = [
    "fetch_gwtc_catalog.py",
    "inspect_catalog_schema.py",
    "build_gwtc_table.py",
    "run_gwtc_log_scan.py",
    "run_gwtc_nulls.py",
    "run_gwtc_bin_stress.py",
    "run_gwtc_variable_stress.py",
    "run_gwtc_threshold_stress.py",
    "run_gwtc_controls.py",
    "make_gwtc_master_table.py",
    "make_gwtc_verdict.py",
]


@pytest.mark.parametrize("name", SCRIPT_NAMES)
def test_script_help_works(name):
    path = os.path.join(SCRIPTS, name)
    assert os.path.exists(path), f"missing script {name}"
    out = subprocess.run([sys.executable, path, "--help"],
                         capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, f"{name} --help failed: {out.stderr}"
    assert "usage" in out.stdout.lower()
