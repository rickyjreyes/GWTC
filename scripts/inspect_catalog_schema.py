#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
inspect_catalog_schema.py
=========================

Inspect the raw GWOSC allevents JSON and report the available fields, their
fill rates, and per-catalog event counts. This audits which admissible
scale-like variables are actually present before any analysis is run.

It does NOT run any WCT scan.
"""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter, defaultdict

DEFAULT_RAW = os.path.join("data", "gwtc_allevents_raw.json")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Inspect GWOSC allevents JSON schema and fill rates.")
    ap.add_argument("--raw", default=DEFAULT_RAW, help="Raw allevents JSON path.")
    ap.add_argument("--top", type=int, default=60, help="Max fields to list.")
    args = ap.parse_args(argv)

    with open(args.raw) as fh:
        payload = json.load(fh)
    events = payload.get("events", {})
    print(f"event entries: {len(events)}")

    cats = Counter(v.get("catalog.shortName") for v in events.values())
    print("\ncatalogs:")
    for c, n in sorted(cats.items(), key=lambda kv: (-kv[1], str(kv[0]))):
        print(f"  {c:28s} {n}")

    # Field fill rates (fraction non-null).
    field_present = defaultdict(int)
    for v in events.values():
        for f, val in v.items():
            if val is not None:
                field_present[f] += 1
    n = max(len(events), 1)
    print("\nfields (fill fraction):")
    for f, c in sorted(field_present.items(), key=lambda kv: -kv[1])[: args.top]:
        print(f"  {f:34s} {c/n:5.2f}  ({c})")

    # Admissible scale-like variables of interest.
    print("\nadmissible scale-like variable availability:")
    for f in [
        "chirp_mass_source",
        "total_mass_source",
        "mass_1_source",
        "mass_2_source",
        "luminosity_distance",
        "redshift",
        "final_mass_source",
        "chi_eff",
        "p_astro",
        "far",
    ]:
        c = field_present.get(f, 0)
        print(f"  {f:24s} {c}/{len(events)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
