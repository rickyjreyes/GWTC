# TODO

Status of the WCT-motivated GWTC log-domain residual diagnostic. Distinguish
**descriptive structure**, **predictive holdout evidence**, and **source/null
discrimination**; a null non-pass is not a causal explanation.

## Done

- [x] Official GWOSC Event API fetch + provenance manifest (URL, date, sha256).
- [x] Clean per-event table with source-frame masses, distance, redshift,
      FAR, p_astro, chi_eff, uncertainties, observing run, source class,
      cumulative dedup flag.
- [x] Historical GWTC-4 168-scan reproduction to floating-point tolerance.
- [x] Dependence-aware GWTC-4 outer nulls N0–N3 with complete-vector global
      calibration.
- [x] Preserve the historical nominal `8.0522 sigma` result as a descriptive
      independent-uniform statistic while separating it from empirical global
      significance.
- [x] Original binned Poisson log-periodic scan, bin/variable/threshold stress,
      jitter and negative controls.
- [x] Strict frozen holdout V1: `PASS_FIXED_MODE`, `p = 0.00129987`.
- [x] V3 non-polynomial/unbinned robustness: training-only Gaussian KDE,
      frozen `k`, amplitude and phase; `PASS_ROBUSTNESS_FIXED_MODE`,
      `p = 0.00089991` on the declared holdout.
- [x] V4 structured astrophysical-null model implemented: continuous broken
      power law + truncated Gaussian peak, no periodic term, with predeclared
      selection exponents `gamma = 0.0, 1.5, 2.5`.
- [x] V4 training-only definition frozen and committed at `k = 9.7`, including
      frozen amplitude and phase. Holdout evaluation is not performed by the
      freeze step.
- [x] Tests for structured-null normalization, deterministic fitting/sampling,
      frozen frequency, and deterministic fixed-statistic null generation.

## Immediate next action

- [ ] **Run the one-shot V4 evaluator exactly as frozen**:

  ```bash
  python -m pytest -q --basetemp=.pytest_tmp
  python scripts/evaluate_gwtc_v4_structured_null.py \
    --null-n 10000 \
    --seed 2718281
  ```

  Then hash and commit `tables/gwtc_v4_structured_null_result.csv` without
  changing `k`, `a,b`, selection exponents, structured-population family, or
  decision threshold.

- [ ] Interpret V4 only as **distinguished / not distinguished from the declared
      non-periodic null family**. Do not translate a non-pass into "ordinary
      astrophysics caused the signal" and do not translate a pass directly into
      WCT confirmation.

## After V4

- [ ] **Full hierarchical population/selection challenge**: propagate event
      posteriors and use injection-calibrated selection effects (or an equivalent
      validated detection-efficiency model) while keeping the residual mode
      externally frozen. This is the strongest conventional source-discrimination
      test still missing.
- [ ] **Independent future-catalog replication**: freeze the exact dimensionless
      mode before a genuinely unseen GWTC release and evaluate once, with no
      retraining on the future catalog. This would carry more evidentiary weight
      than further reanalysis of the already-inspected GWTC-5 holdout.
- [ ] Historical **mode-coherence diagnostic**: instead of only aggregating 168
      marginal p-values, test recurrence/coherence of saved `k_best`, `n_obs`,
      and phase across overlapping subsets using complete-catalog outer nulls.
      Freeze the coherence statistic before generating its null distribution.
- [ ] Posterior-sample-level analysis reusing the legacy `gwtc4_*` HDF5 tooling
      as an independent measurement-uncertainty layer.
- [ ] Per-catalog-version diagnostics (GWTC-4.0 vs GWTC-4.1 vs GWTC-5.0) with
      parameters frozen before each held-out comparison where possible.
- [ ] Add `E_rad` / peak-luminosity variables if/when an official catalog table
      exposes them with sufficient provenance.
- [ ] Seed-stability / Monte Carlo depth checks after the scientific model is
      frozen; use them to assess numerical stability, not to optimize significance.

## Legacy diagnostic maintenance

- [ ] Cap or annotate the binned exploratory `k` grid at the per-binning Nyquist
      limit so aliased high-`k` solutions cannot be mistaken for resolved modes.
- [ ] Preserve family-corrected search accounting for any future exploratory
      scan matrix, while keeping those scans secondary to frozen predictions.
