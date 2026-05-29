# TODO

Status of the WCT-motivated GWTC log-domain residual diagnostic. Diagnostic
only; negative results are retained.

## Done

- [x] Official GWOSC Event API fetch + provenance manifest (URL, date, sha256).
- [x] Clean per-event table with source-frame masses, distance, redshift,
      FAR, p_astro, chi_eff, uncertainties, observing run, source class,
      cumulative dedup flag.
- [x] Core library: Poisson deviance, IRLS GLM, baseline fit, log-periodic
      scan, `n_star`, parametric Poisson-bootstrap null, `p_global`.
- [x] Primary scan + standalone null dump.
- [x] Bin-count stress (20/30/40/50).
- [x] Variable stress (M_chirp, M_total, D_L, M_final, redshift; bounded
      |chi_eff| diagnostic).
- [x] Threshold + observing-run stress (p_astro 0.5/0.9, FAR<1/yr, BBH-only,
      O3a/O3b/O4a/O4b).
- [x] Negative controls (smooth_resample, uniform_ell, jitter, bounded_chi_eff)
      with correct per-type interpretation.
- [x] Master table + a-priori verdict + reliability class.
- [x] Tests: pipeline correctness, overclaim guard, FAIL-tolerance, --help.
- [x] Smoke run on ln(M_chirp); current verdict PARTIAL / Class III.

## Open / future work

- [ ] **More flexible baseline**: the current degree-3 Poisson baseline cannot
      absorb the known astrophysical BBH mass-function peaks, which inflates
      `DeltaD`. Add spline / penalized-spline baselines and confirm whether any
      residual mode survives a baseline that already captures the ~10 and
      ~35 M☉ structure. This is the single most important follow-up before any
      stronger statement.
- [ ] Cap the `k`-grid near the per-binning Nyquist frequency, or report the
      Nyquist line, so aliased high-`k` peaks are not scanned as if physical.
- [ ] Family-corrected `p_family` across the full search matrix (currently the
      search family is recorded but only per-branch `p_global` is reported).
- [ ] Per-catalog-version runs (GWTC-4.0 vs GWTC-4.1 vs GWTC-5.0) to test that
      structure does not appear only because of one catalog update.
- [ ] Posterior-sample-level analysis (per-event posterior draws) reusing the
      legacy `gwtc4_*` HDF5 scripts in the repo root, as an independent layer.
- [ ] Add `E_rad` / peak-luminosity variables if/when a catalog table exposes
      them.
- [ ] Seed-stability sweep (verify verdict is not seed-dependent).
