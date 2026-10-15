#!/usr/bin/env python3
"""Inspect an unmodified GWOSC catalog-index snapshot without touching event data.

A positive entry is only a candidate; release chronology, publication time,
official full cumulative table and physical-event identity require separate
verification before the frozen future manifest can be prepared. Never scores.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from pathlib import Path

SOURCE = "https://gwosc.org/api/v2/catalogs?format=json"
RELEASE = re.compile(r"^GWTC-(\d+)\.(\d+)$")


def audit_index(raw: dict, snapshot_sha256: str) -> dict:
    if not isinstance(raw, dict) or not isinstance(raw.get("results"), list):
        raise ValueError("Not a GWOSC paginated catalog-index response")
    if raw.get("next") is not None or raw.get("page_number") != 1 or raw.get("num_pages") != 1:
        raise ValueError("Incomplete/paginated GWOSC catalog index; collect and verify complete index")
    if len(raw["results"]) != raw.get("results_count"):
        raise ValueError("Catalog index count mismatch")
    releases = {}
    for entry in raw["results"]:
        if not isinstance(entry, dict):
            raise ValueError("Malformed catalog item")
        m = RELEASE.fullmatch(str(entry.get("name", "")))
        if not m:
            continue
        name = entry["name"]
        if name in releases:
            raise ValueError("Duplicate numbered catalog in official snapshot")
        if not str(entry.get("detail_url", "")).startswith(
                "https://gwosc.org/api/v2/catalogs/" + name + "?"):
            raise ValueError("Missing/mismatched official detail URL: " + name)
        releases[name] = (int(m.group(1)), int(m.group(2)))
    if "GWTC-5.0" not in releases:
        raise ValueError("Baseline GWTC-5.0 absent; cannot attest chronology")
    future = sorted((version, name) for name, version in releases.items() if version > (5, 0))
    candidate = future[0][1] if future else None
    return {
        "schema": "GWTC_OFFICIAL_CATALOG_INDEX_SCREEN_V1",
        "source_url_to_archive_independently": SOURCE,
        "snapshot_sha256": snapshot_sha256,
        "numbered_catalogs": sorted(releases, key=lambda n: releases[n]),
        "first_post_5_0_candidate": candidate,
        "status": ("CANDIDATE_REQUIRES_EXTERNAL_RELEASE_ATTESTATION" if candidate
                   else "EXTERNALLY_BLOCKED_NO_POST_5_0_IN_INDEX_SNAPSHOT"),
        "index_is_not_publication_chronology_proof": True,
        "manifest_prepared": False,
        "event_data_read": False,
        "signal_scored": False,
        "next_action": ("Verify official publication chronology/complete cumulative table, "
                        "deduplicated new physical events, timestamp the prepare manifest before scoring. "
                        "Do not treat an API entry or downloaded snapshot alone as authorization.")
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--snapshot", type=Path, required=True,
                   help="Exact official GWOSC catalog-index JSON downloaded without editing")
    p.add_argument("--out", type=Path, required=True,
                   help="New non-scoring read-only evidence file")
    args = p.parse_args(argv)
    if not args.snapshot.is_file():
        p.error("Snapshot missing; download the official catalog-index JSON first")
    if args.out.resolve() == args.snapshot.resolve() or args.out.exists():
        p.error("Refusing to overwrite input or an existing audit record")
    payload = args.snapshot.read_bytes()
    report = audit_index(json.loads(payload), hashlib.sha256(payload).hexdigest())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
