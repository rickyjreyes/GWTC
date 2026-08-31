#!/usr/bin/env python3
"""Calibrate the historical GWTC-4 168-scan population statistic.

This script separates two questions that must not be conflated:

1. Analytic / nominal significance of the observed p-value distribution.
2. Empirical catalog-level significance after preserving dependence among the
   168 overlapping scan statistics.

The second question requires a null matrix produced by rerunning the *entire*
168-scan workflow on synthetic/null catalogs. Each row must be one complete
null catalog and each scan column must correspond to the same scan definition
and ordering as the observed vector.

The historical histogram uses 11 bins with edges
    0, .05, .10, .20, .30, ..., 1.00
so the first two bins have expected count 0.05*N and the remaining bins have
expected count 0.10*N. For N=168 this reproduces expectations 8.4 and 16.8.

This utility does not generate the 168 scan statistics itself. It calibrates
the aggregate statistic once the observed vector and full null-catalog matrix
are available. That separation is deliberate: a matrix made from independent
Uniform(0,1) draws would destroy selector/scan correlations and is not a valid
end-to-end global-null calibration.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import binom, chi2, norm, poisson

HISTORICAL_EDGES = np.array(
    [0.0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00],
    dtype=float,
)
OBSERVED_P_CANDIDATES = ("p_scanmax", "global_p", "p_value", "p")


def validate_pvalues(values: np.ndarray, label: str, expected_n: int | None = None) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim != 1:
        raise SystemExit(f"{label} must be one-dimensional")
    if not np.all(np.isfinite(values)):
        raise SystemExit(f"{label} contains non-finite values")
    if np.any((values < 0.0) | (values > 1.0)):
        raise SystemExit(f"{label} contains values outside [0, 1]")
    if expected_n is not None and values.size != expected_n:
        raise SystemExit(f"{label} contains {values.size} values; expected {expected_n}")
    return values


def load_observed(path: Path, p_column: str | None, expected_n: int) -> tuple[np.ndarray, str]:
    df = pd.read_csv(path)
    if p_column is not None:
        if p_column not in df.columns:
            raise SystemExit(f"Observed p-value column {p_column!r} not found")
        col = p_column
    else:
        matches = [name for name in OBSERVED_P_CANDIDATES if name in df.columns]
        if len(matches) == 1:
            col = matches[0]
        elif len(matches) > 1:
            raise SystemExit(
                "Observed file contains multiple recognized p-value columns; "
                "select one with --p-column"
            )
        elif df.shape[1] == 1:
            col = str(df.columns[0])
        else:
            raise SystemExit(
                "Could not identify observed p-value column. Use --p-column. "
                f"Recognized names: {', '.join(OBSERVED_P_CANDIDATES)}"
            )
    values = validate_pvalues(df[col].to_numpy(float), "Observed p-values", expected_n)
    return values, col


def _candidate_null_columns(df: pd.DataFrame, prefix: str | None) -> list[str]:
    if prefix:
        cols = [str(c) for c in df.columns if str(c).startswith(prefix)]
        if not cols:
            raise SystemExit(f"No null-matrix columns start with prefix {prefix!r}")
        return cols

    cols: list[str] = []
    for c in df.columns:
        s = pd.to_numeric(df[c], errors="coerce")
        if s.notna().all() and ((s >= 0.0) & (s <= 1.0)).all():
            cols.append(str(c))
    return cols


def load_null_matrix(path: Path, expected_n: int, prefix: str | None) -> tuple[np.ndarray, list[str]]:
    df = pd.read_csv(path)
    cols = _candidate_null_columns(df, prefix)
    if len(cols) != expected_n:
        raise SystemExit(
            f"Null matrix has {len(cols)} candidate p-value columns; expected {expected_n}. "
            "Use --null-prefix to select the scan columns explicitly."
        )
    matrix = df[cols].to_numpy(float)
    if matrix.ndim != 2 or matrix.shape[0] == 0:
        raise SystemExit("Null matrix must contain at least one null catalog row")
    if not np.all(np.isfinite(matrix)):
        raise SystemExit("Null matrix contains non-finite values")
    if np.any((matrix < 0.0) | (matrix > 1.0)):
        raise SystemExit("Null matrix contains values outside [0, 1]")
    return matrix, cols


def histogram_chi2(pvalues: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
    counts, _ = np.histogram(pvalues, bins=HISTORICAL_EDGES)
    widths = np.diff(HISTORICAL_EDGES)
    expected = pvalues.size * widths
    statistic = float(np.sum((counts - expected) ** 2 / expected))
    return statistic, counts.astype(int), expected


def safe_z_from_one_sided_p(p: float) -> float:
    if p <= 0.0:
        return math.inf
    if p >= 1.0:
        return -math.inf
    return float(norm.isf(p))


def empirical_upper_tail(null_stats: np.ndarray, observed: float) -> tuple[int, float]:
    ge = int(np.sum(null_stats >= observed))
    p = float((1 + ge) / (1 + null_stats.size))
    return ge, p


def count_below(pvalues: np.ndarray, threshold: float) -> int:
    return int(np.sum(pvalues < threshold))


def aggregate_summary(pvalues: np.ndarray) -> dict[str, object]:
    n = int(pvalues.size)
    chi2_stat, counts, expected = histogram_chi2(pvalues)
    chi2_p = float(chi2.sf(chi2_stat, df=len(HISTORICAL_EDGES) - 2))
    n05 = count_below(pvalues, 0.05)
    n10 = count_below(pvalues, 0.10)

    # Historical paper used Poisson-tail summaries for these threshold counts.
    pois05 = float(poisson.sf(n05 - 1, n * 0.05))
    pois10 = float(poisson.sf(n10 - 1, n * 0.10))

    # Exact marginal count tails are also reported for transparency. They are
    # not a substitute for catalog-level calibration when scans are correlated.
    binom05 = float(binom.sf(n05 - 1, n, 0.05))
    binom10 = float(binom.sf(n10 - 1, n, 0.10))

    return {
        "n_scans": n,
        "histogram_edges": HISTORICAL_EDGES.tolist(),
        "histogram_counts": counts.tolist(),
        "histogram_expected": expected.tolist(),
        "chi2_uniformity": chi2_stat,
        "chi2_df": int(len(HISTORICAL_EDGES) - 2),
        "chi2_analytic_p": chi2_p,
        "chi2_nominal_one_sided_z": safe_z_from_one_sided_p(chi2_p),
        "count_p_lt_0p05": n05,
        "count_p_lt_0p10": n10,
        "poisson_tail_p_lt_0p05": pois05,
        "poisson_tail_p_lt_0p10": pois10,
        "poisson_tail_p_lt_0p05_nominal_z": safe_z_from_one_sided_p(pois05),
        "poisson_tail_p_lt_0p10_nominal_z": safe_z_from_one_sided_p(pois10),
        "binomial_tail_p_lt_0p05": binom05,
        "binomial_tail_p_lt_0p10": binom10,
    }


def null_stat_vectors(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    chi_stats = np.empty(matrix.shape[0], dtype=float)
    n05 = np.empty(matrix.shape[0], dtype=int)
    n10 = np.empty(matrix.shape[0], dtype=int)
    for i, row in enumerate(matrix):
        chi_stats[i] = histogram_chi2(row)[0]
        n05[i] = count_below(row, 0.05)
        n10[i] = count_below(row, 0.10)
    return chi_stats, n05, n10


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Calibrate the historical GWTC-4 168-scan population statistic."
    )
    parser.add_argument("--observed", required=True, help="CSV containing observed scan p-values")
    parser.add_argument("--p-column", default=None, help="Observed p-value column name")
    parser.add_argument(
        "--null-matrix",
        default=None,
        help="CSV: one null catalog per row, one scan p-value per column",
    )
    parser.add_argument(
        "--null-prefix",
        default=None,
        help="Optional prefix selecting p-value columns in the null matrix, e.g. scan_",
    )
    parser.add_argument("--expected-n", type=int, default=168)
    parser.add_argument(
        "--output", default="tables/gwtc4_population_global_result.json"
    )
    args = parser.parse_args()

    if args.expected_n <= 1:
        raise SystemExit("--expected-n must exceed 1")

    observed, observed_column = load_observed(
        Path(args.observed), args.p_column, args.expected_n
    )
    result: dict[str, object] = {
        "schema": "GWTC4_POPULATION_GLOBAL_NULL_V1",
        "observed_file": str(args.observed),
        "observed_p_column": observed_column,
        "historical_statistic_definition": (
            "11-bin Pearson chi-square with edges "
            "[0,.05,.10,.20,.30,.40,.50,.60,.70,.80,.90,1.0]"
        ),
        "observed": aggregate_summary(observed),
        "global_calibration_status": "NOT_CALIBRATED_NO_NULL_MATRIX",
        "interpretation": (
            "Analytic p-to-Z values are nominal until the complete correlated "
            "168-scan workflow is calibrated using null catalogs."
        ),
    }

    if args.null_matrix is not None:
        matrix, columns = load_null_matrix(
            Path(args.null_matrix), args.expected_n, args.null_prefix
        )
        obs = result["observed"]
        assert isinstance(obs, dict)
        chi_null, n05_null, n10_null = null_stat_vectors(matrix)

        chi_ge, chi_emp = empirical_upper_tail(
            chi_null, float(obs["chi2_uniformity"])
        )
        n05_ge, n05_emp = empirical_upper_tail(
            n05_null.astype(float), float(obs["count_p_lt_0p05"])
        )
        n10_ge, n10_emp = empirical_upper_tail(
            n10_null.astype(float), float(obs["count_p_lt_0p10"])
        )

        result["null_matrix_file"] = str(args.null_matrix)
        result["null_catalogs"] = int(matrix.shape[0])
        result["null_scan_columns"] = columns
        result["empirical_catalog_level"] = {
            "chi2_null_ge_observed": chi_ge,
            "chi2_global_p": chi_emp,
            "chi2_global_one_sided_z": safe_z_from_one_sided_p(chi_emp),
            "count_0p05_null_ge_observed": n05_ge,
            "count_0p05_global_p": n05_emp,
            "count_0p05_global_one_sided_z": safe_z_from_one_sided_p(n05_emp),
            "count_0p10_null_ge_observed": n10_ge,
            "count_0p10_global_p": n10_emp,
            "count_0p10_global_one_sided_z": safe_z_from_one_sided_p(n10_emp),
            "monte_carlo_resolution_floor": float(1.0 / (1.0 + matrix.shape[0])),
        }
        result["global_calibration_status"] = "EMPIRICALLY_CALIBRATED_FROM_NULL_MATRIX"
        result["interpretation"] = (
            "Empirical p-values compare the observed aggregate statistics with "
            "complete null-catalog 168-scan vectors, preserving whatever scan "
            "dependence is present in the supplied end-to-end null matrix."
        )

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    obs = result["observed"]
    assert isinstance(obs, dict)
    print(f"Observed scans          : {obs['n_scans']}")
    print(f"Histogram chi2          : {float(obs['chi2_uniformity']):.8g}")
    print(f"Analytic uniformity p   : {float(obs['chi2_analytic_p']):.8g}")
    print(f"Nominal one-sided Z     : {float(obs['chi2_nominal_one_sided_z']):.6g} sigma")
    print(f"Count p < 0.05          : {obs['count_p_lt_0p05']}")
    print(f"Count p < 0.10          : {obs['count_p_lt_0p10']}")
    print(f"Poisson tail p<0.10     : {float(obs['poisson_tail_p_lt_0p10']):.8g}")

    if "empirical_catalog_level" in result:
        emp = result["empirical_catalog_level"]
        assert isinstance(emp, dict)
        print(f"Null catalogs           : {result['null_catalogs']}")
        print(f"Empirical chi2 global p : {float(emp['chi2_global_p']):.8g}")
        print(f"Empirical chi2 global Z : {float(emp['chi2_global_one_sided_z']):.6g} sigma")
        print(f"Empirical count<.10 p   : {float(emp['count_0p10_global_p']):.8g}")
        print(f"MC resolution floor     : {float(emp['monte_carlo_resolution_floor']):.8g}")
    else:
        print("Global calibration      : NOT CALIBRATED (no null matrix supplied)")

    print(f"Wrote                    : {out}")


if __name__ == "__main__":
    main()
