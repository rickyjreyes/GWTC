# GWTC-4 Population Null Sensitivity Study

## Question

The recovered historical GWTC-4 analysis has a real descriptive feature: the
original 168 scan p-values are extremely nonuniform under an independent-
uniform reference (`chi^2 = 95.630952`, nominal `Z = 8.0522 sigma`). The first
correlation-aware catalog null, which fully permutes final-spin medians across
event labels, reduces the complete-vector population significance to about
`1.26 sigma` with 1000 outer catalogs.

Those two facts answer different questions. The purpose of this sensitivity
study is to determine **how much of the nominal 8.05-sigma structure is explained
by overlapping scan geometry, and how much depends on the physical assumptions
built into the full event-label exchangeability null**.

This study does **not** change the observed statistic, scan range, subsets,
thresholds, or historical 168-scan manifest in order to obtain a larger sigma.
Only the predeclared catalog-level null model changes.

## Frozen observed analysis

All null models must use the same recovered historical analysis:

- 86 events;
- 168 scans;
- 56 subset selectors;
- 3 final-spin pairs;
- `k = 0.5 .. 80`, 4000 grid points;
- polynomial baseline degree 2;
- minimum subset size 12;
- top-N selectors `12,16,20,24,32`;
- the saved 168 scan-max `Delta-chi2` values as the observed statistics.

Every sensitivity run first recomputes all historical scan maxima and aborts if
reproduction fails. The established reproduction errors are approximately
`7.1e-15` in `k_best` and `7.4e-09` in `Delta-chi2`.

## Null family

The implementation is
`scripts/generate_gwtc4_population_null_sensitivity.py`.

### N0 — full event-label permutation

`--null-model label`

Finite `final_spin_median` values are freely permuted among all finite-spin
events. Event-specific measurement weights and all non-final-spin quantities
remain fixed.

This is the established exchangeability null. The 1000-catalog result already
found `p_global = 0.1038961`, or `Z = 1.25966 sigma`.

### N1 — mass-stratified permutation

`--null-model mass_stratified`

Final-spin medians are permuted only within predeclared mass quantile strata.
The default is four quantile strata in `log_total_mass_source`.

This preserves broad mass-spin association while breaking fine event-level
ordering. It directly tests whether the full exchangeability null was too
aggressive because it allowed spins to move across very different mass scales.

### N2 — mass + precision stratified permutation

`--null-model mass_precision_stratified`

Final-spin medians are permuted within joint cells formed by mass quantiles and
final-spin uncertainty quantiles. Defaults are four mass bins and two precision
bins.

The historical event-specific weights remain fixed, but donors and recipients
have broadly similar mass and measurement precision. This reduces the concern
that the full event-label null creates unrealistic spin/precision pairings.

### N3 — smooth mass-spin residual permutation

`--null-model mass_residual`

A weighted polynomial baseline is fitted to
`atanh(final_spin) ~ log_total_mass_source` (degree 2 by default). The fitted
residuals are permuted across events and added back to the smooth baseline;
`tanh` maps generated values back into the physical interval `(-1,1)`.

This null explicitly preserves a smooth broad mass-spin relation and asks
whether the historical 168-scan population excess requires fine residual
ordering beyond that smooth relation.

It is still a statistical sensitivity null, not a full LVK astrophysical
population simulation.

## Prespecified first-pass runs

Use 1000 catalogs for each model. The observed statistic remains frozen.
Different seeds are used only to make the runs independently reproducible.

```bash
pytest -q tests/test_population_outer_null.py \
  tests/test_population_null_sensitivity.py \
  tests/test_population_global_null.py

# N1: preserve broad mass association
python scripts/generate_gwtc4_population_null_sensitivity.py \
  --null-model mass_stratified \
  --mass-var log_total_mass_source \
  --mass-bins 4 \
  --outer-n 1000 \
  --seed 20260902

python scripts/run_gwtc_population_global_null.py \
  --observed tables/gwtc4_population_mass_stratified_1000_observed_outercal.csv \
  --p-column outer_global_p \
  --null-matrix tables/gwtc4_population_mass_stratified_1000_null_matrix.csv \
  --output tables/gwtc4_population_mass_stratified_1000_global_result.json

# N2: preserve broad mass and measurement-precision association
python scripts/generate_gwtc4_population_null_sensitivity.py \
  --null-model mass_precision_stratified \
  --mass-var log_total_mass_source \
  --mass-bins 4 \
  --precision-bins 2 \
  --outer-n 1000 \
  --seed 20260903

python scripts/run_gwtc_population_global_null.py \
  --observed tables/gwtc4_population_mass_precision_stratified_1000_observed_outercal.csv \
  --p-column outer_global_p \
  --null-matrix tables/gwtc4_population_mass_precision_stratified_1000_null_matrix.csv \
  --output tables/gwtc4_population_mass_precision_stratified_1000_global_result.json

# N3: preserve a smooth broad mass-spin relation
python scripts/generate_gwtc4_population_null_sensitivity.py \
  --null-model mass_residual \
  --mass-var log_total_mass_source \
  --residual-degree 2 \
  --outer-n 1000 \
  --seed 20260904

python scripts/run_gwtc_population_global_null.py \
  --observed tables/gwtc4_population_mass_residual_1000_observed_outercal.csv \
  --p-column outer_global_p \
  --null-matrix tables/gwtc4_population_mass_residual_1000_null_matrix.csv \
  --output tables/gwtc4_population_mass_residual_1000_global_result.json
```

## Interpretation rules fixed before seeing results

The goal is not to choose whichever null returns the largest sigma. Report all
predeclared nulls together.

| outcome across defensible nulls | interpretation |
|---|---|
| all remain around `p ~ 0.05-0.15` | strong evidence that the apparent 8.05-sigma aggregate is primarily explained by dependence/selector geometry under a broad null family |
| structure-preserving nulls move consistently below `p < 0.01` | evidence that full event-label exchangeability was too destructive; escalate to more physical population nulls |
| one or more predeclared nulls reach `p < 0.0027` (~3 sigma) and replicate across mass coordinates/seeds | meaningful residual population anomaly requiring astrophysical and measurement-aware simulation |
| a null reaches a small p only after changing bins, cuts, variables, or model choices after seeing results | exploratory only; not a discovery significance |
| globally calibrated `p <= 2.87e-7` under a prespecified physically defensible null | statistical 5-sigma threshold, still requiring independent physical interpretation/replication |

A 1000-catalog run cannot directly resolve a 5-sigma tail; its minimum plus-one
p-value is `1/1001`. If a predeclared alternative null produces zero or very few
exceedances, the correct next step is to increase that **same frozen null** or
use tail-calibration methods, not to alter the observed analysis.

## What this study can and cannot establish

If N1-N3 agree with N0, then the original 8.05-sigma number remains an important
descriptive indicator of strong structure in the scan ensemble, but the
structure is not globally rare once the historical dependence architecture is
modeled across several null assumptions.

If N1-N3 materially increase global significance, then the current 1.26-sigma
headline is too specific to the unrestricted event-label exchangeability null.
That would not automatically restore an 8-sigma physical claim; it would show
that a more realistic null family is required before the origin of the base-run
nonuniformity can be considered resolved.

The next level after this sensitivity study is a generative astrophysical null
including population structure, selection effects, and parameter-estimation
uncertainty. That should be built because it is physically justified, not
because it increases sigma.
