#!/usr/bin/env python3
"""Non-scoring future-release readiness gate (#18). Does not inspect residuals.

This cannot certify an official release or someone's historical lack of access.
Such attestations require official provenance and a separately timestamped record.
The frozen future_gwtc_v1.json and one-shot evaluator are NOT modified here.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_future_gwtc import selected_names, check_historical, sha


def check_candidate(table: Path, cfg: dict, release_id: str, source_url: str,
                    published_utc: str) -> dict:
    if not release_id.strip() or release_id.strip().upper() in {
        "GWTC-5", "GWTC-5.0", "GWTC-4", "GWTC-4.0", "GWTC-4.1"
    }:
        raise ValueError("A historical/unspecified release cannot be a future replication")
    if not source_url.startswith("https://") or not published_utc.strip():
        raise ValueError("Official source URL and release timestamp are required")
    if cfg["schema"] != "GWTC_FUTURE_FROZEN_PROTOCOL_V1":
        raise ValueError("Wrong frozen protocol")
    check_historical(cfg)
    future = pd.read_csv(table)
    historical = pd.read_csv(cfg["historical_table_path"])
    selected = selected_names(future, historical, cfg)
    return {
        "schema": "GWTC_FUTURE_NONSCORING_READINESS_V1",
        "status": "CANDIDATE_PENDING_EXTERNAL_FIRST_RELEASE_ATTESTATION",
        "release_id_claimed": release_id,
        "source_url_claimed": source_url,
        "published_utc_claimed": published_utc,
        "table_sha256": sha(table),
        "historical_table_sha256": sha(cfg["historical_table_path"]),
        "candidate_primary_eligible_names": selected.commonName.astype(str).tolist(),
        "candidate_n": int(len(selected)),
        "scores_computed": False,
        "holdout_evaluated": False,
        "required_external_gate": (
            "Independently verify first OFFICIAL cumulative GWTC release after 5.0, "
            "complete source table, publication chronology and new physical identities; "
            "then run prepare, independently timestamp/hash/commit manifest BEFORE one-shot evaluate."
        )
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=Path("protocols/future_gwtc_v1.json"))
    p.add_argument("--candidate-table", type=Path)
    p.add_argument("--release-id")
    p.add_argument("--source-url")
    p.add_argument("--published-utc")
    p.add_argument("--output", type=Path)
    args = p.parse_args(argv)
    if not args.candidate_table:
        payload = {
            "status": "EXTERNALLY_BLOCKED_NO_ATTESTED_UNSEEN_RELEASE",
            "scores_computed": False,
            "protocol": str(args.config),
            "protocol_sha256": sha(args.config),
            "next_action": "Check official chronology; obtain complete first official post-5.0 cumulative table."
        }
    else:
        if not (args.release_id and args.source_url and args.published_utc):
            p.error("Candidate requires --release-id, --source-url, and --published-utc")
        cfg = json.loads(args.config.read_text())
        payload = check_candidate(args.candidate_table, cfg, args.release_id,
                                  args.source_url, args.published_utc)
        payload["protocol_sha256"] = sha(args.config)
    if args.output:
        if args.output.exists():
            p.error("Refusing to overwrite readiness evidence")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
