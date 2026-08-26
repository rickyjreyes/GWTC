#!/usr/bin/env python3
"""Create a frozen chronological GWTC holdout manifest.

The script reads the canonical clean event table and writes an explicit
train/holdout assignment before any holdout residual evaluation is performed.
It prefers catalog/release labels when available and otherwise supports a
user-supplied holdout catalog prefix.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def choose_catalog_column(df: pd.DataFrame) -> str:
    for name in ("catalog_version", "catalog", "catalog_short", "run", "observing_run"):
        if name in df.columns:
            return name
    raise SystemExit(
        "No catalog/run column found. Expected one of: catalog_version, catalog, "
        "catalog_short, run, observing_run."
    )


def choose_event_column(df: pd.DataFrame) -> str:
    for name in ("commonName", "common_name", "event", "event_name", "name"):
        if name in df.columns:
            return name
    raise SystemExit(
        "No event-name column found. Expected one of: commonName, common_name, "
        "event, event_name, name."
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="tables/gwtc_events_clean.csv")
    parser.add_argument("--output", default="tables/gwtc_holdout_manifest.csv")
    parser.add_argument(
        "--holdout-prefix",
        default="GWTC-5",
        help="Catalog/run prefix assigned to holdout (default: GWTC-5).",
    )
    parser.add_argument(
        "--p-astro-min",
        type=float,
        default=0.5,
        help="Apply the same astrophysical-probability floor used by the primary.",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    event_col = choose_event_column(df)
    catalog_col = choose_catalog_column(df)

    if "p_astro" in df.columns:
        df = df[pd.to_numeric(df["p_astro"], errors="coerce") >= args.p_astro_min].copy()

    df[event_col] = df[event_col].astype(str)
    df[catalog_col] = df[catalog_col].astype(str)
    df = df.sort_values([event_col, catalog_col]).drop_duplicates(event_col, keep="last")

    holdout = df[catalog_col].str.startswith(args.holdout_prefix, na=False)
    manifest = pd.DataFrame(
        {
            "event_name": df[event_col],
            "catalog_version": df[catalog_col],
            "split": holdout.map({True: "holdout", False: "train"}),
            "p_astro_min": args.p_astro_min,
            "holdout_rule": f"catalog_prefix:{args.holdout_prefix}",
        }
    ).sort_values(["split", "event_name"])

    if not (manifest["split"] == "holdout").any():
        raise SystemExit(
            f"Holdout rule {args.holdout_prefix!r} selected zero events; refusing to write manifest."
        )
    if not (manifest["split"] == "train").any():
        raise SystemExit("Holdout rule selected all events; refusing to write manifest.")

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(out, index=False)

    counts = manifest["split"].value_counts().to_dict()
    print(f"Wrote {out}: train={counts.get('train', 0)}, holdout={counts.get('holdout', 0)}")


if __name__ == "__main__":
    main()
