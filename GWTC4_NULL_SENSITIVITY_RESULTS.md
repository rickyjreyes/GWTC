# GWTC-4 Prespecified Population-Null Sensitivity Results

This file records the outcomes of the null family frozen in
[`GWTC4_NULL_SENSITIVITY.md`](GWTC4_NULL_SENSITIVITY.md). The observed historical
analysis was not changed between null models: **86 events, 168 scans, 56 subset
selectors, 3 final-spin pairs, `k = 0.5..80` with 4000 grid points, degree-2
baseline, minimum subset size 12, and the same saved scan-max `Delta-chi2`
statistics**.

All three new runs first reproduced the historical scanner with maximum
`|k_recalc-k_saved| = 7.10543e-15` and maximum
`|Delta_recalc-Delta_saved| = 7.43192e-09`; reproduction was **PASS** in every
case. The test suite also passed **15/15** before the sensitivity runs.

## Results

The historical 168-vector remains extremely nonuniform under the naive
independent-uniform reference (`chi^2 = 95.630952`, nominal `Z = 8.0522 sigma`).
The table below reports the scientifically relevant **complete-vector empirical
catalog-null calibration**, not the invalid independent-uniform conversion of
the rank-calibrated p-vectors.

| null | physical structure preserved | outer catalogs | empirical histogram global p | global Z | empirical count `p<0.10` | result |
|---|---|---:|---:|---:|---:|---|
| N0 full event-label permutation | selector/event overlap; non-spin catalog quantities | 1000 | `0.1038961` | `1.25966 sigma` | `0.0629371` | not globally significant |
| N1 mass-stratified permutation | broad total-mass association | 1000 | `0.3446553` | `0.39979 sigma` | `0.1718282` | not globally significant |
| N2 mass+precision stratified permutation | broad total-mass and spin-precision association | 1000 | `0.4285714` | `0.18001 sigma` | `0.1288711` | not globally significant |
| N3 smooth mass-spin residual permutation | smooth degree-2 broad mass-spin relation | 1000 | `0.0639361` | `1.52255 sigma` | `0.0379620` | not globally significant |

For N3, the secondary count-below-0.10 statistic has one-sided
`Z ~= 1.775 sigma`; it is below 0.05 but far from the prespecified escalation
threshold (`p < 0.0027`, about 3 sigma) and is not an independent discovery
statistic.

The rank-calibrated observed p-vectors still produce large **nominal** analytic
numbers if compared with an independent continuous-uniform reference:

| null | nominal analytic histogram Z of rank-calibrated observed p-vector |
|---|---:|
| N1 mass-stratified | `7.121 sigma` |
| N2 mass+precision stratified | `6.310 sigma` |
| N3 mass-residual | `12.040 sigma` |

Those values are diagnostic only. The complete null vectors show directly that
similarly nonuniform dependent p-vectors are common under each corresponding
catalog null.

## Prespecified interpretation applied

No structure-preserving null moved the primary complete-vector histogram test to
`p < 0.01`; none approached `p < 0.0027` (~3 sigma). Instead, the empirical
primary p-values span approximately `0.064..0.429`, with the original unrestricted
permutation at `0.104`.

This outcome supports the first qualitative branch of the prespecified decision
rule, and in N1/N2 is even less anomalous than the originally anticipated
`p ~ 0.05..0.15` range. The unrestricted event-label permutation was therefore
**not uniquely responsible** for reducing the historical 8.05-sigma nominal
aggregate. Preserving broad mass association, mass plus measurement precision,
or a smooth mass-spin relation does not restore a globally significant
population anomaly.

## Scientific conclusion

The **8.0522-sigma base-run number remains a real descriptive fact about the
historical scan output under an independent-uniform reference**. It should not
be erased or described as meaningless. What the null-sensitivity study shows is
that this descriptive nonuniformity is **not globally rare under four materially
different dependence-aware catalog nulls**.

Accordingly, the current defensible statement is:

> The recovered GWTC-4 168-scan ensemble shows striking nominal nonuniformity
> (`8.0522 sigma` under an independent-uniform reference), but the complete-vector
> global significance is not anomalous under the tested dependence-aware null
> family. Across N0-N3 the primary empirical global significance ranges from
> about `0.18 sigma` to `1.52 sigma`.

This does not establish that the underlying catalog contains no physical
structure. It establishes that the historical **population-level >5-sigma
claim is not supported by this scan-aggregate statistic once dependence is
calibrated across the prespecified null family**.

The next scientifically distinct escalation, if pursued, is a generative
astrophysical population null incorporating source-population structure,
selection effects, and parameter-estimation uncertainty. It should be built for
physical fidelity, not selected to increase sigma.

## Generated result files

The local runs produced:

- `tables/gwtc4_population_mass_stratified_1000_observed_outercal.csv`
- `tables/gwtc4_population_mass_stratified_1000_null_matrix.csv`
- `tables/gwtc4_population_mass_stratified_1000_null_delta_matrix.csv`
- `tables/gwtc4_population_mass_stratified_1000_metadata.json`
- `tables/gwtc4_population_mass_stratified_1000_global_result.json`
- `tables/gwtc4_population_mass_precision_stratified_1000_observed_outercal.csv`
- `tables/gwtc4_population_mass_precision_stratified_1000_null_matrix.csv`
- `tables/gwtc4_population_mass_precision_stratified_1000_null_delta_matrix.csv`
- `tables/gwtc4_population_mass_precision_stratified_1000_metadata.json`
- `tables/gwtc4_population_mass_precision_stratified_1000_global_result.json`
- `tables/gwtc4_population_mass_residual_1000_observed_outercal.csv`
- `tables/gwtc4_population_mass_residual_1000_null_matrix.csv`
- `tables/gwtc4_population_mass_residual_1000_null_delta_matrix.csv`
- `tables/gwtc4_population_mass_residual_1000_metadata.json`
- `tables/gwtc4_population_mass_residual_1000_global_result.json`
