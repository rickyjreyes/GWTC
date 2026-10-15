# Remaining GWTC issues: implementation and scientific boundaries (2026-10-14)

This second-stage branch starts from merged PR #17, commit e9dab5b. These changes are run-ready code, data contracts and *read-only* diagnostics; they do NOT invent missing official PE/injections or a genuinely unseen future release. Previously frozen V1/V3/V4/historical N0–N3 data and p-values are unchanged.

| Issue | New deliverable | What still prevents scientific closure |
|---|---|---|
| #15 | `scripts/annotate_legacy_nyquist.py` and committed read-only overlay CSVs for primary and bin stress. **Historical B20 k=32.6975 is ALIASED; do not call it resolved.** Historical verdict fields retained strictly as archival records. | Review presentation and avoid using unannotated historic alias in figures |
| #10 | Full matched observed k/phase/n/event-membership export paired with complete outer-null details; original phase and membership reproduction checked; frozen 1,000 null minimum enforced. | Restore provenance-verified original event-summary CSV and run full outer null |
| #18 | Complete canonical release disposition ledger is hashed and checked between prepare/evaluate, including exclusions and old-event revisions. | First genuinely unseen official post-5.0 cumulative release and external chronology/manifest attestation |
| #11 | `--require-complete`, per-event paired source-frame posterior-vs-point contribution CSV, weighted posterior median and complete coverage diagnostics. | Official event PE HDF5 and independently audited waveform/source-frame units, full 104-event run |
| #8 | Conditional integration of individually validated fixed-hyperparameter hierarchical ratios under externally frozen H0-only baseline posterior weights; mismatch/precision gates. | Full-parameter PE priors, population-model fitting, injection campaigns/exposures, actual posterior/injection runs, independently calibrated scientific inference |
| #12 | Independent checksum-locked cumulative version table ingestion, exact primary identity/diff, official sourced rename mapping, no fake physical deletion. | Supply verified complete official cumulative snapshots and finalize version-sensitivity interpretation |
| #13 | Official definition/unit/frame independent-review preflight before secondary E_rad/peak-luminosity analysis. | Validated official files, source-frame/units checks and new predeclared dependence-corrected scan/null |
| #3 | Machine-readable original acceptance gate; posterior result SHA and complete coverage check; always keeps PARTIAL/Class III. | #11 official PE execution, interpretation, historical provenance caveat, no retroactive promotion |

## Commands after fetching this branch

Run `python -m pytest -q --basetemp=.pytest_tmp`. For archived Nyquist example:
```bash
python scripts/annotate_legacy_nyquist.py --input tables/gwtc_bin_stress.csv \
  --output results/legacy_nyquist/new_bin_stress_annotation.csv
```
Choose a new filename; inputs and existing result files are never overwritten.

For #3, before PE data are available:
```bash
python scripts/gwtc_baseline_v2_acceptance.py --output results/baseline_v2_gate_pre_PE.json
```
For data-dependent procedures see docs/CATALOG_VERSION_AUDIT.md,
POSTERIOR_UNCERTAINTY.md, HIERARCHICAL_CHALLENGE.md,
OFFICIAL_OBSERVABLES.md, HISTORICAL_MODE_COHERENCE.md and FUTURE_REPLICATION.md.
Never execute prospective #18 evaluation against GWTC-5.0 that already influenced V3/V4.


## Targeted #15/#3/#12 follow-up

The publication generator now fail-closes on inconsistent/aliased V2 bin-edge Nyquist, and `scripts/plot_legacy_nyquist.py` renders a separate explicitly aliased B20 archival figure, with non-overwriting output. The frozen baseline-v2 gate now cross-checks actual training CV and V3 KDE provenance, frozen frequencies, split counts, the V2 B40 Nyquist and V3's historical robustness label; it retains PARTIAL/Class III and reports byte-only LF/CRLF representation differences without altering freeze hashes. The original #11 posterior execution remains missing.

For #12, the original GWTC release-tagged rows have zero identity overlap between 5.0 and either 4.x subset. CI reproduced all three archived fixed-mode scores and decomposed the changes. `gwtc_cumulative_frozen_reference.py` is ready for independently supplied official cumulative snapshots but has not been run on them. The added CI evidence job archives the executed #3/#12/#15 local results and SVGs; no original frozen table or figure is overwritten.
