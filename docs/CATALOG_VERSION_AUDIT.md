# GWTC release-row and frozen-mode sensitivity (#12)

The checked-in table contains individual entries tagged GWTC-4.0, GWTC-4.1,
and GWTC-5.0; these are not automatically three independent cumulative
population snapshots. The audit records exact identity overlap, revised fields,
and **not-listed-in-this-version** instead of calling a missing entry a removed
physical merger. It never uses a later version to retune k, phase or bandwidth.

```bash
python scripts/gwtc_catalog_version_audit.py \
 --out-prefix results/gwtc_catalog_version_audit
```

It evaluates the published V4 mode (k=9.7, coefficients and normalizer fixed)
against the fixed V3 KDE reference on version-specific rows with p_astro>=.5
and finite positive source-frame chirp mass; it runs a matching fixed-KDE null
for each sample size. These are **reference-baseline** p-values only, not
independent future-catalog evidence, full population/selection calibration, or
a substitute for an LVK cumulative event manifest with exact canonical
selection/FAR and posterior revisions. Reconcile those differences before
asserting the original #12 scientific gate completed. Preserve all outputs.

## Separate official cumulative release reconciliation

The archived tagged-row audit is not three official cumulative snapshots. Supply three independent, official cumulative CSVs with SHA, official release URL/date and all canonical primary rows via `manifests/cumulative_releases.template.json`. Run:

```bash
python scripts/gwtc_cumulative_manifest_audit.py --manifest local/cumulative_releases.json --out-prefix results/official_cumulative_audit
```

Explicit renames require a source-backed `declared_renames` record; unmatched labels are marked **not listed**, never presumed physically removed. The output is descriptive and hash-locked. The old frozen-V4 vs fixed-V3-KDE reference p-values remain in `results/gwtc_catalog_version_audit_results.json`; they are neither physical population-calibrated nor prospective. Reconcile complete official cumulative field/waveform definitions before #12 scientific closure.


## Exact score decomposition on the archived tagged rows

Run `python scripts/gwtc_version_score_decomposition.py --out-prefix results/local_tagged_decomposition` with a NEW prefix. It verifies the original clean-table and frozen-V4 hashes and independently reconstructs all three published fixed-mode scores without tuning. Its complete physical-identity ledger decomposes each A-to-B statistic change into contributions from **shared IDs**, **IDs not listed in A**, and **IDs not listed in B**, including selection-cut flips on shared IDs. The three terms add exactly to `T_B-T_A`. This is a descriptive comparison of the existing tagged GWOSC rows, NOT three certified official cumulative snapshots, a test of actual physical removal, a posterior/injection-corrected population analysis or a prospective replication. Reference p-values are retained, but differences between p-values are not additive. A verified official cumulative provenance manifest and official parameter/waveform version reconciliation are still necessary for #12 original scientific closure.


## Reproduced tagged-row decomposition (GitHub CI, archival only)

`scripts/gwtc_version_score_decomposition.py` independently reproduced the three archived fixed-mode scores and made the following exact additive comparisons:

| Tagged rows | Common identities | Δ(2logLR) B−A | Shared-ID change | Newly listed in B | Not listed in B |
|---|---:|---:|---:|---:|---:|
| 4.0 → 4.1 | 129 | −0.73806876 | −0.73806876 | 0 | 0 |
| 4.0 → 5.0 | 0 | +7.78995185 | 0 | +9.95018558 | −2.16023373 |
| 4.1 → 5.0 | 0 | +8.52802061 | 0 | +9.95018558 | −1.42216497 |

There are 129, 140 and 161 respective **tagged** entries in the archived clean table. All 129 4.0 IDs recur in 4.1. The 5.0 tagged rows have zero identity overlap with 4.0/4.1. The 11 extra 4.1 rows make zero additional contribution to this fixed score under the stated p_astro/mass cuts. The **p-value differences cannot be decomposed this way** and do not establish prospective replication. This is not proof of physical removal or three independent cumulative catalogs.

### When independently certified cumulative tables are available

After `gwtc_cumulative_manifest_audit.py` validates full official source and canonical identities, run `python scripts/gwtc_cumulative_frozen_reference.py --manifest local/cumulative_releases.json --output results/official_cumulative_fixed_reference_NEW.json`. It rejects output overwrite and invalid frozen V3/V4 provenance; scores each snapshot with exactly the V4 k=9.7 coefficients and the existing V3 training-only KDE reference, using separate predeclared null seeds. It reports canonical-primary counts, exclusions, fixed 2logLR, null exceedances and **reference-KDE p**, explicitly marked neither population/selection corrected nor prospective. Without real official cumulative tables, this runner cannot complete issue #12.

### Frozen Windows text-file provenance

Archived CSV/JSON source SHA256 values reflect historical byte representation. Git's LF checkout may have different raw SHA256 from the same CRLF source. `scripts/gwtc_source_provenance.py` accepts only an EXACT raw SHA or an explicitly reported LF↔CRLF byte transformation matching the historical hash, and records both hashes. Any content or line-order change fails closed. Do not replace or silently regenerate frozen source hashes.
