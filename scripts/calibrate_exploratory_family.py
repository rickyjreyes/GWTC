#!/usr/bin/env python3
"""Dependence-aware calibration of a DECLARED exploratory family (#16).

Input observed CSV: cell_id,statistic. Input null CSV: one complete catalog
realization per row, with a column for EVERY observed cell_id. Simulate the
full selection/scan pipeline upstream; isolated cell bootstraps cannot be
combined by this tool. All statistics are oriented so larger means more extreme.

Use empirical per-column upper tails, then the smallest local tail across the
family; calibrate this against entire dependent null vectors, not independent
p-values. Null-row ranks include that row; this gives an exchangeable
leave-one-out local p based on N-1 other outer catalogs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def calibrate_family(observed: pd.DataFrame, null: pd.DataFrame) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    if list(observed.columns) != ["cell_id", "statistic"]:
        raise ValueError("Observed CSV must have exactly cell_id,statistic columns")
    names = observed.cell_id.astype(str).tolist()
    if not names or len(names) != len(set(names)):
        raise ValueError("Family must be nonempty with unique cell IDs")
    if set(names) != set(null.columns) or len(names) != len(null.columns):
        raise ValueError("Null matrix must include EXACTLY the declared complete family")
    raw = observed.statistic.to_numpy(float)
    matrix = null[names].to_numpy(float)
    if matrix.ndim != 2 or matrix.shape[0] < 2 or not np.isfinite(raw).all() or not np.isfinite(matrix).all():
        raise ValueError("Need at least two complete finite outer-null catalogs and finite observed statistics")
    n = matrix.shape[0]
    # Upper-tail local p with add-one convention for the independent observed catalog.
    p_obs = (1 + (matrix >= raw[None, :]).sum(axis=0)) / (n + 1)
    # For null row i, its own entry is excluded: (1 + count(other >= row_i))/n.
    # Descending maximum rank counts ties conservatively and equals 1+count(other>=row_i).
    p_null = np.column_stack([rankdata(-matrix[:, j], method="max") / n for j in range(len(names))])
    outer_min_p = p_null.min(axis=1)
    p_best = float(p_obs.min())
    family_p = float((1 + (outer_min_p <= p_best).sum()) / (n + 1))
    rows = pd.DataFrame({
        "cell_id": names, "statistic": raw, "local_p": p_obs,
        "family_global_p": family_p, "outer_null_n": n,
        "status": "EXPLORATORY_GLOBAL_CALIBRATION",
    })
    summary = {
        "family_definition": names, "n_cells": len(names),
        "outer_null_n": n, "best_local_p": p_best,
        "family_global_p": family_p,
        "best_cell_ids": [names[i] for i in range(len(names)) if p_obs[i] == p_best],
        "method": "complete-vector empirical min-local-p; descending max rank with ties",
        "evidence_class": "EXPLORATORY_NOT_FROZEN_PREDICTION",
    }
    null_out = pd.DataFrame({"replicate": np.arange(n), "min_local_p": outer_min_p})
    return rows, summary, null_out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--observed", required=True, type=Path)
    p.add_argument("--null-matrix", required=True, type=Path)
    p.add_argument("--out-prefix", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path,
                   help="Predeclared family/provenance JSON; recorded verbatim by SHA256")
    args = p.parse_args()
    import hashlib
    manifest_bytes = args.manifest.read_bytes()
    manifest = json.loads(manifest_bytes)
    if manifest.get("analysis_class") != "exploratory" or not manifest.get("frozen_before_null"):
        raise SystemExit("Manifest must explicitly declare exploratory family frozen before null generation")
    observed = pd.read_csv(args.observed)
    if manifest.get("cell_ids") != observed.cell_id.astype(str).tolist():
        raise SystemExit("Manifest cell_ids must exactly match observed matrix order")
    rows, summary, null_out = calibrate_family(observed, pd.read_csv(args.null_matrix))
    summary["manifest_sha256"] = hashlib.sha256(manifest_bytes).hexdigest()
    summary["observed_sha256"] = hashlib.sha256(args.observed.read_bytes()).hexdigest()
    summary["null_matrix_sha256"] = hashlib.sha256(args.null_matrix.read_bytes()).hexdigest()
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    rows.to_csv(str(args.out_prefix) + "_cells.csv", index=False)
    null_out.to_csv(str(args.out_prefix) + "_outer_null.csv", index=False)
    Path(str(args.out_prefix) + "_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
