#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
GWTC-4 HDF5 Posterior Path Finder
---------------------------------

Purpose
-------
Your previous posterior-only scanner correctly rejected prior paths, but then
found no posterior columns. This utility searches the GWTC-4 HDF5 files deeply
and reports where the real posterior/sample arrays live.

It does NOT run a WCT scan. It only audits HDF5 structure.

It prints and saves:
    outputs_gwtc4_pathfinder/
        gwtc4_hdf5_pathfinder_inventory.csv
        gwtc4_hdf5_pathfinder_hits.csv
        gwtc4_hdf5_pathfinder_summary.json
        gwtc4_hdf5_candidate_prefixes.csv

Run
---
    python gwtc4_find_posterior_paths.py --input "*combined_PEDataRelease.hdf5" --limit-files 3

Then inspect:
    outputs_gwtc4_pathfinder/gwtc4_hdf5_pathfinder_hits.csv

If posterior paths exist, you should see names containing:
    posterior
    posterior_samples
    samples/posterior
    publication
    posterior_samples/mass_1_source
    posterior_samples/chirp_mass_source

If only priors/samples appears, then these files do not contain posterior samples
in normal PE-release form, or the posterior samples are stored in compound/table
datasets that need special handling.
"""

from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path
from typing import Any, Dict, List, Tuple

import h5py
import numpy as np
import pandas as pd


OUT_DIR_DEFAULT = "outputs_gwtc4_pathfinder"

PARAM_KEYWORDS = [
    "mass_1", "mass1", "m1",
    "mass_2", "mass2", "m2",
    "chirp", "mc",
    "total_mass", "mtotal",
    "mass_ratio", "q",
    "symmetric_mass_ratio", "eta",
    "chi_eff", "chieff",
    "chi_p", "chip",
    "luminosity", "distance",
    "redshift",
    "posterior",
    "sample",
    "samples",
    "publication",
]

BAD_KEYWORDS = [
    "prior",
    "priors",
    "config_file",
    "calibration",
    "psd",
    "injection",
    "skymap",
]

GOOD_PREFIX_HINTS = [
    "posterior",
    "posterior_samples",
    "samples/posterior",
    "publication",
    "pesummary/posterior",
    "samples",
]

COMMON_PARAM_NAMES = [
    "mass_1_source",
    "mass_2_source",
    "mass_1",
    "mass_2",
    "chirp_mass_source",
    "chirp_mass",
    "total_mass_source",
    "total_mass",
    "mass_ratio",
    "symmetric_mass_ratio",
    "chi_eff",
    "chi_p",
    "luminosity_distance",
    "redshift",
]


def safe_name(path: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", path).strip("_")


def discover_files(patterns: List[str], limit_files: int | None = None) -> List[str]:
    files = []
    for pat in patterns:
        hits = glob.glob(pat)
        if hits:
            files.extend(hits)
        elif Path(pat).exists():
            files.append(pat)
    files = sorted(set(files))
    if limit_files is not None and limit_files > 0:
        files = files[:limit_files]
    return files


def is_numeric_dtype(dtype: Any) -> bool:
    try:
        return bool(np.issubdtype(dtype, np.number))
    except Exception:
        return False


def is_compound_dtype(dtype: Any) -> bool:
    try:
        return dtype.fields is not None
    except Exception:
        return False


def dtype_field_names(dtype: Any) -> List[str]:
    if not is_compound_dtype(dtype):
        return []
    return list(dtype.fields.keys())


def norm_path(name: str) -> str:
    return name.lower().replace("\\", "/")


def classify_path(name: str, dtype: str = "", fields: str = "") -> Dict[str, Any]:
    p = norm_path(name)
    d = str(dtype).lower()
    f = str(fields).lower()

    has_param_keyword = any(k in p or k in f for k in PARAM_KEYWORDS)
    has_good_hint = any(k in p for k in GOOD_PREFIX_HINTS)
    has_bad_hint = any(k in p for k in BAD_KEYWORDS)
    has_posterior = "posterior" in p
    has_prior = "prior" in p or "priors" in p
    has_sample = "sample" in p or "samples" in p

    # This is not a final policy; it is a ranking score for inspection.
    score = 0
    if has_posterior:
        score += 100
    if has_good_hint:
        score += 20
    if has_sample:
        score += 10
    if has_param_keyword:
        score += 10
    if has_bad_hint:
        score -= 100
    if has_prior:
        score -= 100

    return {
        "has_param_keyword": bool(has_param_keyword),
        "has_good_hint": bool(has_good_hint),
        "has_bad_hint": bool(has_bad_hint),
        "has_posterior": bool(has_posterior),
        "has_prior": bool(has_prior),
        "has_sample": bool(has_sample),
        "candidate_score": int(score),
    }


def inventory_file(path: str, max_items: int | None = None) -> pd.DataFrame:
    rows = []
    path = str(path)

    with h5py.File(path, "r") as f:
        def visitor(name: str, obj: Any):
            if max_items is not None and len(rows) >= max_items:
                return

            kind = "group" if isinstance(obj, h5py.Group) else "dataset" if isinstance(obj, h5py.Dataset) else type(obj).__name__

            shape = ""
            dtype = ""
            n_rows = np.nan
            n_cols = np.nan
            fields = []

            if isinstance(obj, h5py.Dataset):
                shape_obj = getattr(obj, "shape", None)
                dtype_obj = getattr(obj, "dtype", None)
                shape = str(shape_obj)
                dtype = str(dtype_obj)
                if shape_obj is not None and len(shape_obj) >= 1:
                    try:
                        n_rows = int(shape_obj[0])
                    except Exception:
                        n_rows = np.nan
                if shape_obj is not None and len(shape_obj) >= 2:
                    try:
                        n_cols = int(shape_obj[1])
                    except Exception:
                        n_cols = np.nan
                if dtype_obj is not None:
                    fields = dtype_field_names(dtype_obj)

            cls = classify_path(name, dtype=dtype, fields=" ".join(fields))

            rows.append({
                "file": path,
                "name": name,
                "kind": kind,
                "shape": shape,
                "dtype": dtype,
                "n_rows": n_rows,
                "n_cols": n_cols,
                "fields": ",".join(fields),
                "is_numeric": bool(isinstance(obj, h5py.Dataset) and is_numeric_dtype(getattr(obj, "dtype", None))),
                "is_compound": bool(isinstance(obj, h5py.Dataset) and is_compound_dtype(getattr(obj, "dtype", None))),
                **cls,
            })

        f.visititems(visitor)

    return pd.DataFrame(rows)


def read_preview(path: str, dataset: str, max_rows: int = 5) -> Dict[str, Any]:
    """
    Try to preview a dataset without loading huge arrays.
    Handles numeric, string, and compound/table datasets.
    """
    out: Dict[str, Any] = {
        "file": path,
        "dataset": dataset,
        "preview_ok": False,
        "preview_error": "",
        "preview": "",
        "columns": "",
    }

    try:
        with h5py.File(path, "r") as f:
            d = f[dataset]
            dtype = d.dtype
            shape = d.shape

            if dtype.fields is not None:
                fields = list(dtype.fields.keys())
                out["columns"] = ",".join(fields)
                arr = d[: min(max_rows, shape[0])]
                preview_rows = []
                for row in arr:
                    preview_rows.append({field: scalar_to_jsonable(row[field]) for field in fields[:30]})
                out["preview"] = json.dumps(preview_rows, default=str)[:2000]
            else:
                arr = d[: min(max_rows, shape[0])] if len(shape) else d[()]
                out["preview"] = np.array2string(np.asarray(arr), threshold=50)[:2000]

            out["preview_ok"] = True
    except Exception as exc:
        out["preview_error"] = str(exc)

    return out


def scalar_to_jsonable(x: Any) -> Any:
    try:
        if isinstance(x, bytes):
            return x.decode(errors="replace")
        if np.isscalar(x):
            y = x.item()
            if isinstance(y, bytes):
                return y.decode(errors="replace")
            if isinstance(y, (float, int, str, bool)) or y is None:
                return y
            return str(y)
        return str(x)
    except Exception:
        return str(x)


def candidate_prefixes(hits: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize likely group prefixes that contain parameter datasets.
    """
    rows = []
    if hits.empty:
        return pd.DataFrame()

    for _, r in hits.iterrows():
        name = str(r["name"])
        parts = name.split("/")
        for depth in range(1, len(parts)):
            prefix = "/".join(parts[:depth])
            rows.append({
                "file": r["file"],
                "prefix": prefix,
                "name": name,
                "candidate_score": r.get("candidate_score", 0),
                "has_posterior": r.get("has_posterior", False),
                "has_prior": r.get("has_prior", False),
            })

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    agg = df.groupby(["file", "prefix"], as_index=False).agg(
        n_hits=("name", "nunique"),
        max_score=("candidate_score", "max"),
        posterior_hits=("has_posterior", "sum"),
        prior_hits=("has_prior", "sum"),
    )
    agg = agg.sort_values(["max_score", "posterior_hits", "n_hits"], ascending=[False, False, False])
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", nargs="*", default=["*combined_PEDataRelease.hdf5"], help="HDF5 glob(s)")
    ap.add_argument("--out-dir", default=OUT_DIR_DEFAULT)
    ap.add_argument("--limit-files", type=int, default=3, help="Limit files inspected; use 0 for all")
    ap.add_argument("--max-items", type=int, default=None, help="Optional max HDF5 objects per file")
    ap.add_argument("--preview", action="store_true", help="Preview top hit datasets")
    ap.add_argument("--preview-n", type=int, default=20, help="Number of top hit datasets to preview")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    limit = None if args.limit_files == 0 else args.limit_files
    files = discover_files(args.input, limit_files=limit)
    if not files:
        raise FileNotFoundError(f"No HDF5 files found for {args.input}")

    print(f"[files] inspecting {len(files)} file(s)")
    for f in files:
        print(" ", f)

    inventories = []
    for i, f in enumerate(files, 1):
        print(f"[inventory] {i}/{len(files)} {f}")
        inv = inventory_file(f, max_items=args.max_items)
        inventories.append(inv)

    inventory = pd.concat(inventories, ignore_index=True)
    inv_path = out_dir / "gwtc4_hdf5_pathfinder_inventory.csv"
    inventory.to_csv(inv_path, index=False)

    hits = inventory[
        (inventory["kind"] == "dataset")
        & (
            inventory["has_param_keyword"]
            | inventory["has_posterior"]
            | inventory["has_sample"]
            | inventory["is_compound"]
        )
    ].copy()

    hits = hits.sort_values(
        ["candidate_score", "has_posterior", "has_prior", "n_rows"],
        ascending=[False, False, True, False],
        na_position="last",
    ).reset_index(drop=True)

    hits_path = out_dir / "gwtc4_hdf5_pathfinder_hits.csv"
    hits.to_csv(hits_path, index=False)

    prefixes = candidate_prefixes(hits)
    prefixes_path = out_dir / "gwtc4_hdf5_candidate_prefixes.csv"
    prefixes.to_csv(prefixes_path, index=False)

    previews = []
    if args.preview and not hits.empty:
        for _, r in hits.head(args.preview_n).iterrows():
            previews.append(read_preview(str(r["file"]), str(r["name"])))
        preview_df = pd.DataFrame(previews)
        preview_df.to_csv(out_dir / "gwtc4_hdf5_pathfinder_previews.csv", index=False)

    summary = {
        "input": args.input,
        "files_inspected": files,
        "n_inventory_rows": int(len(inventory)),
        "n_hit_rows": int(len(hits)),
        "n_candidate_prefix_rows": int(len(prefixes)),
        "top_hits": hits.head(50).to_dict(orient="records"),
        "top_prefixes": prefixes.head(50).to_dict(orient="records") if not prefixes.empty else [],
        "verdict": {
            "posterior_hits": int(hits["has_posterior"].sum()) if not hits.empty else 0,
            "prior_hits": int(hits["has_prior"].sum()) if not hits.empty else 0,
            "compound_hits": int(hits["is_compound"].sum()) if not hits.empty else 0,
            "numeric_hits": int(hits["is_numeric"].sum()) if not hits.empty else 0,
        },
    }

    with open(out_dir / "gwtc4_hdf5_pathfinder_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print("\n[summary]")
    print(json.dumps(summary["verdict"], indent=2))

    print("\n[top hits]")
    show_cols = ["name", "kind", "shape", "dtype", "fields", "candidate_score", "has_posterior", "has_prior", "has_sample"]
    if not hits.empty:
        print(hits.head(40)[show_cols].to_string(index=False))
    else:
        print("No hits found.")

    print("\n[top candidate prefixes]")
    if not prefixes.empty:
        print(prefixes.head(40).to_string(index=False))
    else:
        print("No prefixes found.")

    print("\n[save]")
    print(inv_path)
    print(hits_path)
    print(prefixes_path)
    print(out_dir / "gwtc4_hdf5_pathfinder_summary.json")
    if args.preview:
        print(out_dir / "gwtc4_hdf5_pathfinder_previews.csv")


if __name__ == "__main__":
    main()
