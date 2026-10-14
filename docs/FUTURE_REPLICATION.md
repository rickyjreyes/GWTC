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

The prepare phase now records a complete non-residual `release_ledger`: primary old/new identities, p_astro exclusions, invalid masses, and historical IDs not listed as primary. The evaluate phase verifies the exact ledger and historical table hash in addition to the existing frozen config and table hashes. **An external auditor must still verify official cumulative completeness, identity revisions, release chronology and pre-inspection timestamp.** Do not execute #18 on GWTC-5.0 or an old-event-only revision.
