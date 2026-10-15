#!/usr/bin/env python3
"""One-shot, no-refit future GWTC protocol (#9).

Prepare a metadata-only manifest first and commit/hash it externally. Evaluate
only afterward using --locked-manifest-sha256, with immutable scientific config.
The future release MUST be the first official cumulative GWTC release after 5.0.
This condition requires an auditable release attestation; code cannot prove
that an operator has never looked at the results outside this process.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent))
from gwtc_unbinned_kde import fixed_log_likelihood_ratio, fixed_model_null

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def primary_mask(df):
    if "is_primary_entry" not in df:
        raise ValueError("Missing official canonical-primary indicator")
    return df.is_primary_entry.astype(str).str.strip().str.lower().isin(("true", "1", "yes"))

def selected_names(future, historical, cfg):
    required = {"commonName", "M_chirp", "p_astro", "is_primary_entry"}
    if not required.issubset(future):
        raise ValueError("Missing required columns: " + str(sorted(required - set(future))))
    if future.commonName.isna().any() or historical.commonName.isna().any():
        raise ValueError("Missing event name")
    primary = future.loc[primary_mask(future)].copy()
    if primary.commonName.duplicated().any():
        raise ValueError("Duplicate primary commonName; fix provenance before locking")
    prior = set(historical.commonName.astype(str))
    new = primary.loc[~primary.commonName.astype(str).isin(prior)].copy()
    new = new.loc[pd.to_numeric(new.p_astro, errors="coerce") >= cfg["selection"]["p_astro_min_inclusive"]]
    if not len(new):
        raise ValueError("No genuinely unseen events meeting frozen p_astro cut")
    return new.sort_values("commonName")

def release_ledger(future, historical, cfg):
    """Metadata-only dispositions for EVERY official canonical primary identity.

    An external auditor must still attest official completeness/release order.
    No residual, frequency or statistic is computed in this function.
    """
    required={"commonName","M_chirp","p_astro","is_primary_entry"}
    if not required.issubset(future):raise ValueError("Incomplete official cumulative table")
    if future.commonName.isna().any() or historical.commonName.isna().any():
        raise ValueError("Missing merger identity")
    primary=future.loc[primary_mask(future)].copy()
    if primary.commonName.duplicated().any():raise ValueError("Duplicate canonical primary merger")
    old=set(historical.commonName.astype(str))
    new=primary.loc[~primary.commonName.astype(str).isin(old)].copy()
    prob=pd.to_numeric(new.p_astro,errors="coerce")
    mask=prob>=cfg["selection"]["p_astro_min_inclusive"]
    eligible=new.loc[mask]
    mass=pd.to_numeric(eligible.M_chirp,errors="coerce").to_numpy(float)
    return {
      "n_release_rows":int(len(future)),
      "n_primary_rows":int(len(primary)),
      "old_primary_names":sorted(set(primary.commonName.astype(str))&old),
      "new_primary_names_before_cut":sorted(new.commonName.astype(str).tolist()),
      "new_names_below_or_missing_pastro":sorted(new.loc[~mask,"commonName"].astype(str).tolist()),
      "eligible_names":sorted(eligible.commonName.astype(str).tolist()),
      "eligible_invalid_mass_names":sorted(eligible.loc[~(np.isfinite(mass)&(mass>0)),"commonName"].astype(str).tolist()),
      "historical_names_not_listed_as_primary":sorted(old-set(primary.commonName.astype(str)))}

def check_historical(cfg):
    path = Path(cfg["historical_table_path"])
    if sha(path) != cfg["historical_table_sha256"]:
        raise ValueError("Historical catalog hash changed")
    if sha(cfg["v3_model_path"]) != cfg["v3_model_sha256"]:
        raise ValueError("Frozen training KDE changed")
    v4 = json.loads(Path(cfg["v4_mode_path"]).read_text())
    for name, expected in cfg["frozen_mode"].items():
        got = v4["fixed_mode"][name]
        if not np.isclose(got, expected, rtol=0, atol=1e-12):
            raise ValueError("V4 fixed mode differs from frozen future config: " + name)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage", choices=("prepare", "evaluate"))
    p.add_argument("--future-table", required=True, type=Path,
                   help="Complete official cumulative table; never a hand-picked subset")
    p.add_argument("--config", default="protocols/future_gwtc_v1.json", type=Path)
    p.add_argument("--source-url", help="Official source URL (prepare only)")
    p.add_argument("--release-id", help="Official first post-5.0 release name (prepare only)")
    p.add_argument("--published-utc", help="Official date of release (prepare only)")
    p.add_argument("--manifest", required=True, type=Path)
    p.add_argument("--locked-manifest-sha256", help="Previously archived/committed manifest SHA256 (evaluate only)")
    p.add_argument("--output", default="tables/gwtc_future_replication_result.json", type=Path)
    args = p.parse_args()
    cfg = json.loads(args.config.read_text())
    if cfg["schema"] != "GWTC_FUTURE_FROZEN_PROTOCOL_V1": raise SystemExit("Wrong protocol")
    check_historical(cfg)
    future = pd.read_csv(args.future_table)
    historical = pd.read_csv(cfg["historical_table_path"])
    new = selected_names(future, historical, cfg)
    ledger = release_ledger(future, historical, cfg)
    if args.stage == "prepare":
        if not (args.source_url and args.release_id and args.published_utc):
            raise SystemExit("Official source URL, release ID and publication date required")
        if args.manifest.exists(): raise SystemExit("Refusing to overwrite an existing manifest")
        obj = {
            "schema": "GWTC_FUTURE_EVENT_MANIFEST_V1",
            "assertion": "first official cumulative GWTC release after GWTC-5.0; verify externally",
            "release_id": args.release_id, "source_url": args.source_url,
            "published_utc": args.published_utc, "config_sha256": sha(args.config),
            "table_sha256": sha(args.future_table),
            "historical_table_sha256": sha(cfg["historical_table_path"]),
            "selected_event_names": new.commonName.astype(str).tolist(),
            "release_ledger": ledger,
            "eligible_rows_before_finite_mass_filter": int(len(new)),
        }
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        args.manifest.write_text(json.dumps(obj, sort_keys=True, indent=2) + "\n")
        print("Commit/archive manifest BEFORE any residual evaluation; SHA256:", sha(args.manifest))
        return
    if not args.locked_manifest_sha256 or sha(args.manifest) != args.locked_manifest_sha256:
        raise SystemExit("Locked manifest SHA256 mismatch or missing")
    m = json.loads(args.manifest.read_text())
    if (m["schema"] != "GWTC_FUTURE_EVENT_MANIFEST_V1"
        or m["config_sha256"] != sha(args.config)
        or m["table_sha256"] != sha(args.future_table)
        or m["selected_event_names"] != new.commonName.astype(str).tolist()
        or m.get("release_ledger") != ledger
        or m.get("historical_table_sha256") != sha(cfg["historical_table_path"])):
        raise SystemExit("Frozen manifest does not match current release/config/event selection")
    mass = pd.to_numeric(new.M_chirp, errors="coerce").to_numpy(float)
    valid = np.isfinite(mass) & (mass > 0)
    ell = np.log(mass[valid])
    if not len(ell): raise SystemExit("No valid source-frame positive masses")
    model = json.loads(Path(cfg["v3_model_path"]).read_text())
    f = cfg["frozen_mode"]; n = cfg["null"]["n"]; seed = cfg["null"]["seed"]
    observed = fixed_log_likelihood_ratio(ell, f["k"], f["a_cos"], f["b_sin"], f["normalizer_Z"])
    null = fixed_model_null(kde_centers=model["train_ell"], bandwidth=float(model["selected_bandwidth"]),
        n_events=len(ell), k=f["k"], a=f["a_cos"], b=f["b_sin"], normalizer=f["normalizer_Z"],
        null_n=n, seed=seed)
    exceed = int(np.sum(null >= observed)); pval = float((1 + exceed)/(1 + n))
    result = {
       "schema":"GWTC_FUTURE_ONE_SHOT_RESULT_V1", "status":"PROSPECTIVE_IF_MANIFEST_PRECOMMITTED",
       "config_sha256":sha(args.config), "manifest_sha256":sha(args.manifest),
       "table_sha256":sha(args.future_table), "release_id":m["release_id"],
       "selected_events":m["selected_event_names"], "n_eligible":len(new),
       "n_valid_mass":len(ell), "n_invalid_mass":int((~valid).sum()),
       "n_outside_training_span":int(((ell < model["train_ell_min"])|(ell > model["train_ell_max"])).sum()),
       "statistic_2loglr":observed,"null_n":n,"seed":seed,"exceedances":exceed,"p":pval,
       "verdict":"PASS_FUTURE_FROZEN_MODE" if observed>0 and pval<cfg["null"]["threshold"] else "FAIL_FUTURE_FROZEN_MODE",
       "interpretation":"Frozen KDE reference, not full posterior/injection-corrected population inference"
    }
    if args.output.exists(): raise SystemExit("Refusing to overwrite one-shot result")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
