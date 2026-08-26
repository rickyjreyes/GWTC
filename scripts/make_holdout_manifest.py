#!/usr/bin/env python3
"""Create a frozen chronological GWTC holdout manifest.

By default this uses the canonical deduplicated ``is_primary_entry`` rows from
``gwtc_events_clean.csv`` before applying the p_astro cut.  This keeps the
holdout population identical to the cumulative primary analysis and avoids
selecting a different catalog copy of the same physical event.
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


def build_manifest(
    df: pd.DataFrame,
    holdout_prefix: str = "GWTC-5",
    p_astro_min: float = 0.5,
    primary_only: bool = True,
) -> pd.DataFrame:
    event_col = choose_event_column(df)
    catalog_col = choose_catalog_column(df)

    work = df.copy()
    if primary_only and "is_primary_entry" in work.columns:
        primary = work["is_primary_entry"]
        if primary.dtype != bool:
            primary = primary.astype(str).str.lower().isin(("true", "1", "yes"))
        work = work[primary].copy()

    if "p_astro" in work.columns:
        work = work[pd.to_numeric(work["p_astro"], errors="coerce") >= p_astro_min].copy()

    work[event_col] = work[event_col].astype(str)
    work[catalog_col] = work[catalog_col].astype(str)

    # Canonical primary rows should already be unique.  Keep this guard for
    # synthetic/legacy tables and fail if ambiguity survives the canonical cut.
    if work[event_col].duplicated().any():
        duplicated = work.loc[work[event_col].duplicated(keep=False), event_col].unique()
        raise ValueError(
            "Holdout input contains duplicate event names after canonical filtering: "
            + ", ".join(map(str, duplicated[:5]))
        )

    holdout = work[catalog_col].str.startswith(holdout_prefix, na=False)
    manifest = pd.DataFrame(
        {
            "event_name": work[event_col],
            "catalog_version": work[catalog_col],
            "split": holdout.map({True: "holdout", False: "train"}),
            "p_astro_min": p_astro_min,
            "holdout_rule": f"catalog_prefix:{holdout_prefix}",
            "primary_only": primary_only,
        }
    ).sort_values(["split", "event_name"])

    if not (manifest["split"] == "holdout").any():
        raise ValueError(
            f"Holdout rule {holdout_prefix!r} selected zero events; refusing to write manifest."
        )
    if not (manifest["split"] == "train").any():
        raise ValueError("Holdout rule selected all events; refusing to write manifest.")
    return manifest


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
    parser.add_argument(
        "--all-entries",
        action="store_true",
        help="Do not restrict to is_primary_entry. Not recommended for the strict holdout.",
    )
    args = parser.parse_args()

    df = pd.read_csv(args.input)
    try:
        manifest = build_manifest(
            df,
            holdout_prefix=args.holdout_prefix,
            p_astro_min=args.p_astro_min,
            primary_only=not args.all_entries,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(out, index=False)

    counts = manifest["split"].value_counts().to_dict()
    print(
        f"Wrote {out}: train={counts.get('train', 0)}, "
        f"holdout={counts.get('holdout', 0)}, canonical_primary={not args.all_entries}"
    )


if __name__ == "__main__":
    main()
