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
      frozen amplitude and phase.
- [x] V4 one-shot evaluation completed with `Delta 2logL = 9.95018558` and
      **`PASS_ALL_STRUCTURED_NULLS`**:
      - `gamma = 0.0`: `55/10000`, `p = 0.0055994401`
      - `gamma = 1.5`: `1/10000`, `p = 0.0001999800`
      - `gamma = 2.5`: `2/10000`, `p = 0.0002999700`
      - worst-case `p = 0.0055994401`
- [x] V4 result preserved in `tables/gwtc_v4_structured_null_result.csv`.
- [x] Tests for structured-null normalization, deterministic fitting/sampling,
      frozen frequency, phase-specific signal definition, and deterministic
      fixed-statistic null generation.

## Current interpretation

The V4 result shows that the exact frozen `k = 9.7`, amplitude, and phase
statistic is uncommon under **all three declared non-periodic structured
population scenarios**. This strengthens the fixed-mode GWTC evidence beyond
the original binned scan and beyond V3's KDE robustness test.

It does **not** replace a full LVK hierarchical population/selection analysis,
and GWTC-5 had already been inspected in earlier V2/V3 work, so V4 is a
robustness/source-discrimination challenge rather than an independent future
catalog replication.

## Immediate next action

- [ ] **Full hierarchical population/selection challenge**: propagate event
      posteriors and use injection-calibrated selection effects (or an equivalent
      validated detection-efficiency model) while keeping the residual mode
      externally frozen. This is now the strongest conventional
      source-discrimination test still missing.
- [ ] **Independent future-catalog replication**: freeze the exact dimensionless
      mode before a genuinely unseen GWTC release and evaluate once, with no
      retraining on the future catalog. This would carry more evidentiary weight
      than further reanalysis of the already-inspected GWTC-5 holdout.

## Additional follow-up

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
