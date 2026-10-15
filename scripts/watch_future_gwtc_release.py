#!/usr/bin/env python3
"""Archive a fresh official GWOSC index snapshot without touching event data.

A numbered post-5.0 entry triggers only an external-review alert; neither this
program nor its scheduled workflow downloads event tables, prepares a holdout
manifest or invokes evaluate_future_gwtc.py. #18 remains open until the
original one-shot scientific test is legitimately executed.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
from urllib.request import Request, urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_official_future_catalog_index import SOURCE, audit_index


def screen_snapshot(payload: bytes, *, observed_utc: str) -> dict:
    raw = json.loads(payload)
    report = audit_index(raw, hashlib.sha256(payload).hexdigest())
    report["observed_utc_not_publication_time"] = observed_utc
    report["workflow_action"] = (
        "EXTERNAL_CHRONOLOGY_AND_EVENT_IDENTITY_REVIEW_REQUIRED"
        if report["first_post_5_0_candidate"] else "WAIT_FOR_OFFICIAL_POST_5_0_INDEX_ENTRY"
    )
    report["first_release_attested"] = False
    report["future_manifest_sha256"] = None
    report["future_result"] = None
    return report


def write_snapshot(out_dir: Path, payload: bytes, *, observed_utc: str) -> dict:
    if out_dir.exists():
        raise FileExistsError("Refusing to overwrite a prior catalog-index snapshot directory")
    # Parse and validate fully before persisting either artifact.
    report = screen_snapshot(payload, observed_utc=observed_utc)
    out_dir.mkdir(parents=True)
    (out_dir / "official_gwosc_index_raw.json").write_bytes(payload)
    (out_dir / "non_scoring_screen.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return report


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", required=True, type=Path, help="Fresh immutable-by-convention artifact directory")
    p.add_argument("--snapshot", type=Path, help="Offline test/review of exact previously downloaded bytes")
    args = p.parse_args(argv)
    if args.out_dir.exists():
        p.error("Output directory exists; never replace earlier observation evidence")
    if args.snapshot:
        payload = args.snapshot.read_bytes()
    else:
        req = Request(SOURCE, headers={"User-Agent": "GWTC-future-release-nonscoring-watch/1.0",
                                        "Accept": "application/json"})
        with urlopen(req, timeout=30) as response:
            final_url = response.geturl()
            if final_url != SOURCE:
                raise ValueError("Unexpected redirect or index source URL: " + final_url)
            payload = response.read(5 * 1024 * 1024 + 1)
            if len(payload) > 5 * 1024 * 1024:
                raise ValueError("Index exceeds fixed 5 MiB limit")
    stamp = datetime.now(timezone.utc).isoformat()
    report = write_snapshot(args.out_dir, payload, observed_utc=stamp)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
