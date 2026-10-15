# Historical binned Nyquist presentation audit (#15)

**Maintenance annotation only.** PR #17 already enforces a common observed/null
resolved grid in new exploratory scans. Historical binned tables and figures
remain unchanged. A previously published maximum above its own binning's
Nyquist frequency is **not** a resolved mode merely because it was selected
or assigned a small conditional p-value.

For uniform log-bin centers, `k_N = pi * B / Delta_ell_active`; the exact
boundary is excluded. The archival `n_star = k_best * Delta_ell_active/(2*pi)`
lets the diagnostic recover the active support of each existing CSV row.
Run from repository root:

```bash
python scripts/annotate_legacy_nyquist.py
# output: results/legacy_nyquist_annotations.csv
pytest -q tests/test_legacy_nyquist_annotations.py
```

It reads `gwtc_scan_primary.csv`, `gwtc_bin_stress.csv`, and
`gwtc_master_results.csv` without modifying their values, original
verdicts, figures, null scores, or commits. `RESOLVED` only means below the
sampling limit; it is **not** independent evidence or global significance.
`INSUFFICIENT_METADATA` means no physical inference is permitted.

## Required caption on reused historical presentation

| Archived source | B | Historical k* | Classification |
| --- | ---: | ---: | --- |
| Primary M_chirp scan | 20 | 32.69747899159664 | **ALIASED_OR_NYQUIST** (Nyquist approximately 14.14) |
| Bin stress M_chirp | 20 | 32.69747899159664 | **ALIASED_OR_NYQUIST** |
| Bin stress M_chirp | 30/40/50 | 4.483193277310924 | Below their respective Nyquist limits |

**Caption to retain wherever the original 20-bin maximum is reused:**
"Legacy 20-bin exploratory maximum k≈32.70 is above the corresponding
log-bin Nyquist bound k_N≈14.14, hence aliased/unresolved; the historical
statistics are preserved for provenance, not promoted as evidence of a resolved
high-frequency feature. Modern scans limit observed and simulated grids
identically below Nyquist. Unbinned V3/V4 tests are distinct."

The original `figures/scan_M_chirp_B30.png` is a legacy display, **not**
an unmodified primary prospective evidence figure; its underlying per-k
curve is not reconstructed from a single saved image. Do not paint an
invented Nyquist crossing onto that image. The publication-facing
`figures/publication/` V2/V3/V4 suite is a separate frozen-mode chain and
must not be mislabeled as the historical exploratory bin scan. The sampling
limit does not apply to unbinned likelihood analyses.

Presentation acceptance: verify generated `results/legacy_nyquist_annotations.csv`
and its alias labels and use this caption/table in any publication reusing
legacy high-k figures. No reclassification of historical science is performed
by this metadata-only audit.
