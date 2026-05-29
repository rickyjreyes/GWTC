#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
build_gwtc_table.py
===================

Build a clean per-event table from the raw GWOSC allevents JSON.

Each output row is one event-entry from one catalog version. We keep
catalog_version as a column so each catalog can be analyzed separately, and we
mark a deduplicated "primary entry" per commonName so a cumulative-confident
view can be selected without double counting events that appear in multiple
catalog releases.

Admissible scale-like variables emitted (source frame where applicable):
    M_chirp     = chirp_mass_source
    M_total     = total_mass_source  (fallback m1_source + m2_source)
    M_final     = final_mass_source
    D_L         = luminosity_distance
    redshift    = redshift
    chi_eff     = chi_eff            (bounded; diagnostic only, not a scale test)

Selection / quality columns: p_astro, far, observing run, crude source class.

Source-class heuristic (declared, not authoritative): using a ~3 Msun neutron
star cutoff on source-frame component masses:
    BBH  : m1 > 3 and m2 > 3
    NSBH : m1 > 3 and m2 <= 3
    BNS  : m1 <= 3 and m2 <= 3
This is a coarse mass-based label, NOT an LVK classification.

Output: tables/gwtc_events_clean.csv
"""

from __future__ import annotations

import argparse
import json
import math
import os

import pandas as pd

DEFAULT_RAW = os.path.join("data", "gwtc_allevents_raw.json")
DEFAULT_OUT = os.path.join("tables", "gwtc_events_clean.csv")

# Approximate observing-run GPS boundaries (LVK standard run windows).
RUN_WINDOWS = [
    ("O1", 1126051217, 1137254417),
    ("O2", 1164556817, 1187733618),
    ("O3a", 1238166018, 1253977218),
    ("O3b", 1256655618, 1269363618),
    ("O4a", 1368921618, 1389466818),
    ("O4b", 1396796418, 1430000000),
]

# Catalog preference for cumulative dedup (higher index = more preferred).
CATALOG_PREFERENCE = [
    "Initial_LIGO_Virgo",
    "O1_O2-Preliminary",
    "GWTC-1-marginal",
    "GWTC-2",
    "GWTC-2.1-marginal",
    "GWTC-2.1-auxiliary",
    "GWTC-3-marginal",
    "O3_IMBH_marginal",
    "GWTC-1-confident",
    "GWTC-2.1-confident",
    "GWTC-3-confident",
    "O3_Discovery_Papers",
    "O4_Discovery_Papers",
    "GWTC-4.0",
    "GWTC-4.1",
    "GWTC-5.0",
]

NS_MAX_MSUN = 3.0


def assign_run(gps):
    if gps is None or (isinstance(gps, float) and math.isnan(gps)):
        return "unknown"
    for name, lo, hi in RUN_WINDOWS:
        if lo <= gps < hi:
            return name
    return "other"


def source_class(m1, m2):
    if m1 is None or m2 is None:
        return "unknown"
    try:
        m1 = float(m1)
        m2 = float(m2)
    except (TypeError, ValueError):
        return "unknown"
    if math.isnan(m1) or math.isnan(m2):
        return "unknown"
    if m1 > NS_MAX_MSUN and m2 > NS_MAX_MSUN:
        return "BBH"
    if m1 > NS_MAX_MSUN and m2 <= NS_MAX_MSUN:
        return "NSBH"
    return "BNS"


def chirp_from_components(m1, m2):
    if m1 is None or m2 is None:
        return None
    try:
        m1 = float(m1)
        m2 = float(m2)
    except (TypeError, ValueError):
        return None
    if m1 <= 0 or m2 <= 0:
        return None
    return (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2


def catalog_rank(cat):
    try:
        return CATALOG_PREFERENCE.index(cat)
    except ValueError:
        return -1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Build clean per-event GWTC table from raw GWOSC JSON.")
    ap.add_argument("--raw", default=DEFAULT_RAW, help="Raw allevents JSON path.")
    ap.add_argument("--out", default=DEFAULT_OUT, help="Clean events CSV output path.")
    args = ap.parse_args(argv)

    with open(args.raw) as fh:
        payload = json.load(fh)
    events = payload.get("events", {})

    rows = []
    for key, v in events.items():
        m1 = v.get("mass_1_source")
        m2 = v.get("mass_2_source")
        mtot = v.get("total_mass_source")
        if mtot is None and m1 is not None and m2 is not None:
            try:
                mtot = float(m1) + float(m2)
            except (TypeError, ValueError):
                mtot = None
        mchirp = v.get("chirp_mass_source")
        if mchirp is None:
            mchirp = chirp_from_components(m1, m2)
        gps = v.get("GPS")
        rows.append({
            "event_key": key,
            "commonName": v.get("commonName"),
            "catalog_version": v.get("catalog.shortName"),
            "GPS": gps,
            "observing_run": assign_run(gps),
            "source_class": source_class(m1, m2),
            "M_chirp": mchirp,
            "M_total": mtot,
            "mass_1_source": m1,
            "mass_2_source": m2,
            "M_final": v.get("final_mass_source"),
            "D_L": v.get("luminosity_distance"),
            "redshift": v.get("redshift"),
            "chi_eff": v.get("chi_eff"),
            "p_astro": v.get("p_astro"),
            "far": v.get("far"),
            # asymmetric uncertainties for jitter controls
            "M_chirp_lower": v.get("chirp_mass_source_lower"),
            "M_chirp_upper": v.get("chirp_mass_source_upper"),
            "M_total_lower": v.get("total_mass_source_lower"),
            "M_total_upper": v.get("total_mass_source_upper"),
            "D_L_lower": v.get("luminosity_distance_lower"),
            "D_L_upper": v.get("luminosity_distance_upper"),
        })

    df = pd.DataFrame(rows)

    # Mark cumulative primary entry: highest-preference catalog per commonName.
    df["_rank"] = df["catalog_version"].map(catalog_rank)
    df["is_primary_entry"] = False
    idx = df.groupby("commonName")["_rank"].idxmax()
    df.loc[idx, "is_primary_entry"] = True
    df = df.drop(columns=["_rank"])

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    df.to_csv(args.out, index=False)

    n_primary = int(df["is_primary_entry"].sum())
    n_primary_pastro = int(
        ((df["is_primary_entry"]) & (df["p_astro"].fillna(0) >= 0.5)).sum()
    )
    print(f"[build] wrote {args.out}: {len(df)} entries, "
          f"{df['commonName'].nunique()} unique events")
    print(f"[build] primary (deduplicated) entries: {n_primary}")
    print(f"[build] primary with p_astro>=0.5: {n_primary_pastro}")
    print("[build] source-class counts (all entries):")
    print(df["source_class"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
