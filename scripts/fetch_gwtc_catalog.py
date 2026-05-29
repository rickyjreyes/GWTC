#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
fetch_gwtc_catalog.py
=====================

Fetch official GWTC catalog event metadata from the public GWOSC Event API.

Source (official LVK / GWOSC):
    https://gwosc.org/eventapi/json/allevents/

This downloads the canonical per-event summary table (source-frame masses,
luminosity distance, redshift, FAR, p_astro, chi_eff, with asymmetric
uncertainties) for every catalog GWOSC publishes, including GWTC-4.x and
GWTC-5.0. We use the official API rather than scraped secondary summaries.

Outputs:
    data/gwtc_allevents_raw.json   (raw API payload, verbatim)
    data/manifest.csv              (provenance: URL, access date, sha256, counts)

Provenance (URL, access date, sha256, event counts) is recorded so the build
is auditable and reproducible.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import sys
import urllib.request
from collections import Counter

DEFAULT_URL = "https://gwosc.org/eventapi/json/allevents/"
DEFAULT_OUT = os.path.join("data", "gwtc_allevents_raw.json")
DEFAULT_MANIFEST = os.path.join("data", "manifest.csv")


def sha256_of_bytes(b: bytes) -> str:
    h = hashlib.sha256()
    h.update(b)
    return h.hexdigest()


def fetch(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "wct-gwtc-diagnostic/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Fetch official GWOSC GWTC event catalog (allevents JSON).",
    )
    ap.add_argument("--url", default=DEFAULT_URL, help="GWOSC event API URL.")
    ap.add_argument("--out", default=DEFAULT_OUT, help="Raw JSON output path.")
    ap.add_argument("--manifest", default=DEFAULT_MANIFEST, help="Provenance manifest CSV.")
    ap.add_argument("--timeout", type=int, default=120, help="HTTP timeout seconds.")
    ap.add_argument(
        "--offline-ok",
        action="store_true",
        help="If set and download fails but --out exists, reuse the cached file.",
    )
    args = ap.parse_args(argv)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    access_date = _dt.datetime.now(_dt.timezone.utc).isoformat()
    raw = None
    try:
        print(f"[fetch] GET {args.url}", file=sys.stderr)
        raw = fetch(args.url, timeout=args.timeout)
    except Exception as exc:  # network/HTTP error
        print(f"[fetch] download failed: {exc}", file=sys.stderr)
        if args.offline_ok and os.path.exists(args.out):
            print(f"[fetch] reusing cached {args.out}", file=sys.stderr)
            with open(args.out, "rb") as fh:
                raw = fh.read()
        else:
            return 2

    # Validate JSON and count.
    try:
        payload = json.loads(raw.decode("utf-8"))
        events = payload.get("events", {})
    except Exception as exc:
        print(f"[fetch] invalid JSON payload: {exc}", file=sys.stderr)
        return 3

    digest = sha256_of_bytes(raw)
    with open(args.out, "wb") as fh:
        fh.write(raw)

    cats = Counter(v.get("catalog.shortName") for v in events.values())
    print(f"[fetch] wrote {args.out} ({len(raw)} bytes, sha256={digest[:16]}...)", file=sys.stderr)
    print(f"[fetch] {len(events)} event entries across {len(cats)} catalogs", file=sys.stderr)

    # Write provenance manifest.
    os.makedirs(os.path.dirname(args.manifest) or ".", exist_ok=True)
    import csv
    with open(args.manifest, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["key", "value"])
        w.writerow(["source", "GWOSC Event API (official LVK / GWOSC)"])
        w.writerow(["url", args.url])
        w.writerow(["access_date_utc", access_date])
        w.writerow(["raw_file", args.out])
        w.writerow(["sha256", digest])
        w.writerow(["n_bytes", len(raw)])
        w.writerow(["n_event_entries", len(events)])
        for cat, n in sorted(cats.items(), key=lambda kv: (-kv[1], str(kv[0]))):
            w.writerow([f"catalog::{cat}", n])
    print(f"[fetch] wrote provenance manifest {args.manifest}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
