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
