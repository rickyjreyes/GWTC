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
