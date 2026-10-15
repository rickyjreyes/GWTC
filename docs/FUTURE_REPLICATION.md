# First truly unseen GWTC release — frozen prospective protocol (#9)

The exact published future scientific rule is in
`protocols/future_gwtc_v1.json`. It freezes the *V4* exact k=9.7 amplitude,
phase and normalizer and the *V3* KDE reference including its original training
centers and bandwidth. It freezes the first official **cumulative GWTC release
after GWTC-5.0**, the new-event-only selection, canonical primary rows,
source-frame chirp mass in solar masses, p_astro>=0.5, no scan or model refit,
10,000 Monte Carlo null replicates and seed 20261014; p<0.05 with positive
statistic is the declared test. All eligible new events with valid positive
source-frame masses are included, including outside the training domain
(reported separately); invalid mass rows are counted and reported.

When an eligible release exists, obtain and archive the complete official
cumulative table, provenance and date. Run:

```bash
python scripts/evaluate_future_gwtc.py prepare \
 --future-table official_cumulative.csv \
 --source-url OFFICIAL_URL --release-id OFFICIAL_ID \
 --published-utc YYYY-MM-DD --manifest results/future_manifest.json
sha256sum results/future_manifest.json
git add results/future_manifest.json && git commit -m "Lock future GWTC identities before residual evaluation"
```

Externally timestamp/archive the manifest commit and SHA before residual
evaluation. Then run exactly once:

```bash
python scripts/evaluate_future_gwtc.py evaluate \
 --future-table official_cumulative.csv --manifest results/future_manifest.json \
 --locked-manifest-sha256 ARCHIVED_SHA256 \
 --output results/future_result.json
```

The evaluator refuses changed source tables, selection lists, config, input
hashes and result overwrite. It cannot prevent someone inspecting results
outside the executable; provenance and precommit timestamps must be audited
externally. A release containing only revised OLD events is not a new holdout.
V1/V3/V4 already saw GWTC-5.0; they do not become independent just because
this protocol was added. A negative future result remains negative.


## Execution readiness and blocking condition (#18; 2026-10-15)

Issue #18 cannot be scientifically closed by this code change: **there is no
independently attested first official cumulative post-GWTC-5.0 release in this
repository and no new-event manifest has been timestamped.** The original
freeze file, event selection, frozen V4 parameters, V3 KDE and once-only
evaluator remain unchanged. As a non-scoring preflight run:

```bash
python scripts/check_future_gwtc_readiness.py --config protocols/future_gwtc_v1.json
# prints EXTERNALLY_BLOCKED_NO_ATTESTED_UNSEEN_RELEASE; no statistic
```

Only after an eligible release has been independently identified, prepare a
candidate **without residual inspection** using the complete official table:

```bash
python scripts/check_future_gwtc_readiness.py \
 --candidate-table official_cumulative.csv \
 --release-id OFFICIAL_NEW_RELEASE --source-url OFFICIAL_URL \
 --published-utc YYYY-MM-DD --output results/future_readiness.json
```

This records input hashes and candidate new primary event identities, not
mass values, likelihood ratios or p-values. The candidate remains
`CANDIDATE_PENDING_EXTERNAL_FIRST_RELEASE_ATTESTATION`; a string supplied
as a release ID is not verification. Independently verify the official release
chronology, complete cumulative official table, event identity/version history,
physical-merger deduplication, p_astro cut and all new eligible event inclusions.
Archive/timestamp this verification. Only THEN follow the original
`prepare` → external timestamp + immutable manifest hash/commit →
single `evaluate --locked-manifest-sha256` workflow above. Neither #18 nor
the future-result verdict may be marked complete from preflight alone.

## Official release index screen (additional non-scoring evidence)

The GWOSC official catalog index is `https://gwosc.org/api/v2/catalogs?format=json`.
At the 2026-10-15 check it lists GWTC-5.0 but no numbered post-5.0
release. Archive the *raw* official JSON without editing and screen it:

```bash
curl -fsSL 'https://gwosc.org/api/v2/catalogs?format=json' -o results/gwosc_catalog_index_raw.json
python scripts/audit_official_future_catalog_index.py \
  --snapshot results/gwosc_catalog_index_raw.json \
  --out results/gwosc_index_screen.json
```

The screen rejects incomplete pagination, count mismatch, missing baseline,
duplicate names and non-GWOSC detail links. It records the input SHA256 and
either `EXTERNALLY_BLOCKED_NO_POST_5_0_IN_INDEX_SNAPSHOT` or
`CANDIDATE_REQUIRES_EXTERNAL_RELEASE_ATTESTATION`. This is an index snapshot
check only, **not proof of chronological first publication, a prepared
manifest, new physical events, or test authorization**. Check official
publication dates and complete event-table provenance independently before
using the existing frozen prepare/evaluate steps. Index screening never reads
event masses or invokes the evaluator. Do not use low-latency alerts or
revisions to old events as a substitute for the frozen official release rule.
