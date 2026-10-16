# Official release reconciliation: issue #12

Baseline: PR #19 branch at `09b3aabd77dc40d98df06dee529dfb420f3591fa`; PR #17 (`49b9f7849c6533605acf835f151cdcddd4ce52c9`) was merged and its numerical/provenance machinery is reused. The execution PR is #23; no open predecessor is assumed merged.

Three official GWOSC catalog JSON snapshots were retrieved on 2026-10-15 and are preserved with hashes. Preferred parameter-estimation metadata were checked on 2026-10-15/16 for **all 278 selected version/event rows**. Every selected value exactly matches the official preferred `chirp_mass_source` in `M_sun`; no detector-frame chirp mass was substituted. Waveform assumptions are recorded per event:

| Release | Raw entries | Selected positive source masses | Preferred PE families |
|---|---:|---:|---|
| 4.0 | 129 | 86 | 55 C00:Mixed; 29 C00:Mixed+XO4a; 1 C00:IMRPhenomXPHM-SpinTaylor; 1 C00:Mixed:HighSpin |
| 4.1 | 140 | 88 | 55 C00:Mixed; 29 C00:Mixed+XO4a; 3 C00:IMRPhenomXPHM-SpinTaylor; 1 C00:Mixed:HighSpin |
| 5.0 | 161 | 104 | 103 C00:IMRPhenomXPHM-SpinTaylor; 1 C01:Mixed |

These are metadata checks, not proof that distinct waveform families give equivalent posteriors. In particular, the [4.0 release documentation](https://gwosc.org/GWTC-4.0/) describes the original likelihood problem, later reweighting and rerun preferred PE. The retrieved 4.0 portal table is therefore an **as-retrieved official snapshot**, not a claim to have recovered the original publication-time bytes. [4.1](https://gwosc.org/GWTC-4.1/) reanalyzes O4a; the [5.0 portal release](https://gwosc.org/GWTC-5.0/) documents O4b and links separately to the current cumulative event list.

## Identity and release scope

All 129 4.0 identities appear in 4.1, with 11 additional listed candidates. Shared names have identical GPS times. There are no duplicate physical-GPS entries within a release and no same-GPS unmatched-name rename candidates across these snapshots. Zero names overlap 4.0/4.1 with 5.0. The disjointness describes different observing-period lists; it is not physical removal of all O4a mergers.

The current official [cumulative JSON endpoint](https://gwosc.org/eventapi/json/GWTC/) was also retrieved: 391 entries, assigned to 5.0 (161), 4.1 (140), 2.1-confident (54), 3-confident (35), and 1-confident (1). Its SHA-256 is `9802865cc8b41508f4073769758f7dbaf5498b22a907849d567083635e7e862a`. This confirms a separate current cumulative product. It does **not** provide three immutable historical all-run snapshots as of the 4.0, 4.1 and 5.0 publication dates. We do not create a union of current rows and label it an official historical snapshot.

## Frozen numerical execution

The prior numerical evaluator and additive decomposition are reused unchanged, with the new output explicitly labeled incremental. Training KDE, k=9.7, phase coefficients, normalizer, p_astro>=0.5, positive-mass selection, 10,000 nulls and seed 20261014 (fixed per-version offsets) are preserved. No outcome-driven retuning occurred.

| Official tagged rows | Selected n | Frozen 2 log LR | Null exceedances / 10,000 | Conditional KDE reference p |
|---|---:|---:|---:|---:|
| GWTC-4.0 | 86 | 2.1602337300 | 494 | 0.0494950505 |
| GWTC-4.1 | 88 | 1.4221649663 | 586 | 0.0586941306 |
| GWTC-5.0 | 104 | 9.9501855796 | 5 | 0.0005999400 |

These reproduce the earlier tagged-row numbers. The 4.0→4.1 score change is -0.7380687638, entirely from shared identities (including two selection flips); the 11 newly listed rows contribute zero under the frozen positive-mass cut. The ledger retains each event's mass, selection and signed score contribution. The 4.x→5.0 comparisons are between disjoint tagged samples, not revision effects on the same mergers.

These nulls are the frozen training-KDE reference experiment. They do not calibrate a hierarchical population/selection test, PE uncertainty, cross-version predictive independence or multiplicity across exploratory release comparisons. The result is retrospective catalog sensitivity; it does not supply a new unseen prediction.

## Criteria and remaining obligation

| Criterion | Verified status |
|---|---|
| Official manifests and event identity/diff | Executed for all three incremental snapshots; current cumulative snapshot separately identified |
| Reproducible frozen comparisons | Executed; original reference numbers reproduced, no parameter refit |
| Exact mass/PE definition | Source-mass units and preferred-waveform metadata verified for all 278 scored rows; waveform-family differences retained |
| Added/renamed/revised/deduplicated entries | Exact name/GPS reconciliation; field revisions and selection flips in ledgers; no silent aliases or removals |
| Historical cumulative version reconciliation | **OPEN:** three authoritative all-run as-of manifests, their version-selection policy and original PE revision provenance are not available in the retrieved inputs |
| Sensitivity versus prediction verdict | Retrospective, conditional reference sensitivity; no promotion of historical PARTIAL/Class III |

Classification: **PARTIAL_OPEN**. A complete original cumulative-revision interpretation needs authoritative historical snapshots; present release lists cannot be relabeled to satisfy it.

## Reproduction and evidence

```bash
python scripts/gwtc_official_version_execution.py --sources results/official_2026-10-15 --output NEW/version_execution
python -m pytest -q tests/test_official_version_execution.py
```

Use a new output directory. Frozen-input checks accept only the existing documented LF/CRLF byte transformation. `results/official_2026-10-15/version_execution/` contains three manifests, `event_diff.csv`, `score_change_ledger.csv`, and `results.json` with source hashes, preferred-PE identities and frozen provenance. The metadata snapshot archive preserves exact downloaded JSON bytes. Tests reject incorrect units, duplicate GPS identities, shared-name identity drift and silent rename merging; missing masses remain in the manifests. The overall execution receipt hashes these artifacts and scripts.
