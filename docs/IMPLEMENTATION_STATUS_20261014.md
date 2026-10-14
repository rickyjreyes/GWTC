# GWTC issue implementation status — 2026-10-14

This branch proposes *code and preregistration*, not newly measured physics.
No existing V1/V3/V4, historical 168-scan, N0–N3, manifest or raw data
artifacts are intentionally changed. Every new data-dependent output requires
a separate execution and provenance review.

| Issue | New implementation | Status before external runs |
|---|---|---|
| #15 | Common per-binning Nyquist filter, identical observed/null grid, reject aliased direct scans, tests, metadata | Code submitted; run CI and inspect plot/table effects |
| #16 | Complete-catalog dependence-preserving full-family rank calibration, fail on missing cells, tests | Code submitted; supply new predeclared exploratory family |
| #9 | Future V4 mode/V3 baseline frozen JSON and two-stage prepare/evaluate, hashed new-event-only manifest | Preparation implemented; unseen release still pending |
| #11 | Official per-event HDF5 SHA/dataset/weights ingestion, posterior one-draw-per-event sensitivity diagnostic | Code submitted; download and verify all appropriate official event products; NOT a calibrated posterior likelihood |
| #12 | Version-specific event diff and fixed-mode reference diagnostic | Code submitted; full cumulative version matching and calibrated comparison still pending |
| #8 | Rate-marginalized fixed-mode importance reweighting core with numerical gates | Core only; externally validate p0(theta), PE priors and injection proposals, compute actual likelihood challenge |
| #14 | Predeclared V1/V3/V4 multi-seed/depth runner and exact-binomial uncertainty | Prepared; long numerical matrix not yet run |
| #13 | Official HDF5 observable presence audit, never substitute unofficial estimates | Metadata gate only; official field/units audit and secondary analysis pending |
| #10 | Frozen historical coherence statistic and optional full k/phase/n/selector-membership export from complete outer-null generator, tests | Recovered historical event summary and 1,000-catalog null replay still pending |
| #3 | Original acceptance checklist reconciled with V1/V3/V4; no changed criteria | Official posterior uncertainty run still pending; preserve Class III |

No GitHub issues should be closed merely because this PR is merged. Preserve
original outcome labels, including negative/null results. Particularly V3 and
V4 are robustness checks on the already-inspected GWTC-5 sample, not independent
prospective tests.

## Data release provenance

Official documentation: https://gwosc.org/GWTC-5.0/
GWTC-5.0 PE: https://zenodo.org/records/20348005 and
https://zenodo.org/records/20348006 .
Cumulative search sensitivity: https://zenodo.org/records/19500052 .
These multi-GB inputs are NOT embedded or fetched in CI; scripts require local
provenance manifests and checksum verification.

## Reproduce the code checks

```bash
python -m pip install -r requirements.txt
python -m compileall -q scripts tests
python -m pytest -q --basetemp=.pytest_tmp
```

Do not execute `evaluate_future_gwtc.py evaluate` until a *future* release
exists and the corresponding prepare manifest has been committed/timestamped.
