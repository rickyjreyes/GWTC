# RESULTS

WCT-motivated GWTC log-domain residual diagnostic. **This is a diagnostic
harness only — it does not prove WCT and does not replace LVK population
inference.** All results below are reproducible from the committed protocols,
model freezes, manifests, and result tables.

## Current evidence hierarchy

The repository contains distinct tests answering different questions. They
should not be collapsed into a single sigma value or interpreted as though a
null model provides a unique physical explanation.

| layer | result | current reading |
|---|---|---|
| Historical GWTC-4 168-scan aggregate | nominal `8.0522 sigma` under independent-uniform p-vector reference | Real descriptive nonuniformity; independence assumption is not a valid catalog-level calibration. |
| GWTC-4 dependence-aware N0–N3 outer nulls | empirical global `0.18–1.52 sigma` on the same p-vector histogram statistic | The **historical aggregate statistic** is not globally calibrated >5 sigma under these nulls. This is not a proof that the catalog structure is nonphysical. |
| Strict frozen holdout V1 | `holdout_p = 0.00129987`, `PASS_FIXED_MODE` | A training-frozen chirp-mass mode predicts the declared holdout under the frozen polynomial-baseline protocol. |
| Unbinned KDE V3 | `Delta 2logL = 10.0122`, `8/10000` null exceedances, `p = 0.00089991`, `PASS_ROBUSTNESS_FIXED_MODE` | The frozen mode survives removal of histogram binning and replacement of the polynomial baseline by a training-only Gaussian KDE. GWTC-5 was already inspected in V2, so this is robustness rather than a new independent catalog replication. |
| Structured astrophysical-null V4 | **frozen; result not yet committed** | Tests the exact frozen `k = 9.7`, amplitude, and phase against non-periodic broken-power-law + Gaussian-peak population models with predeclared selection weighting. |

Committed fixed-mode results:

- `tables/gwtc_frozen_holdout_result.csv`
- `tables/gwtc_v3_unbinned_kde_holdout_result.csv`
- `tables/gwtc_v4_frozen_structured_null.json` (frozen V4 definition; evaluation pending)

The current inferential emphasis is therefore **prediction and exact mode
survival**, not the historical count of small p-values.

A critical interpretive rule for the repository is:

> **Not distinguished from a declared null family** means that the tested
> statistic is not diagnostic between the observed catalog and that null. It
> does **not** prove that the null's physical mechanism caused the observation.

Likewise, a structured-null pass would show that the exact frozen mode is
uncommon under that declared non-periodic family; it would not by itself prove
WCT or replace a full LVK hierarchical population/selection analysis.

## Historical GWTC-4 population-level aggregate result

The restored historical analysis contains **168 scans = 56 subset selectors x
3 final-spin pairs**. The original scan table is committed as
`tables/gwtc4_population_observed.csv`, and the historical scanner is preserved
as `scripts/gwtc4_wct_subset_scan_compound.py`.

### Nominal historical calculation

Using the actual 168 historical `scan_null_p` values:

| population statistic | literal recomputation | nominal one-sided Z |
|---|---:|---:|
| distribution nonuniformity | `chi^2 = 95.630952`, `p = 4.06578e-16` | **8.0522 sigma** |
| count at `p < 0.10` | 47 / 168; Poisson `p = 1.16e-9`; exact Binomial `p = 4.85e-11` | **5.97 sigma Poisson / 6.47 sigma exact-binomial** |
| count at `p < 0.05` | 16 / 168; Poisson `p = 0.0125`; exact Binomial `p = 0.0105` | about **2.24-2.31 sigma** |

The historical paper text describes "11 equal bins," but its figure and expected
counts correspond to `[0,.05]`, `[.05,.10]`, then nine width-`.10` bins. The
audit uses that figure-compatible construction.

These numbers are **nominal marginal/independence-based diagnostics**. The 168
scans reuse the same 86 events across heavily overlapping ranking and quantile
subsets, so treating the full p-vector as independent uniform draws is not an
adequate global null.

### Restored-scan verification

Before generating any population nulls, the recovered implementation was forced
to reproduce all historical scan maxima from the saved event summary:

| verification quantity | result |
|---|---:|
| events | 86 |
| observed scans | 168 |
| unique subsets | 56 |
| scan pairs | 3 |
| dynamic final-spin subsets | 7 |
| max `|k_recalc-k_saved|` | `7.10543e-15` |
| max `|Delta_recalc-Delta_saved|` | `7.43192e-09` |
| reproduction | **PASS** |

The historical selector/scanner architecture is therefore reproducible from the
restored files to floating-point tolerance.

### N0 — unrestricted event-label permutation

`scripts/generate_gwtc4_population_null_matrix.py` generates whole-catalog outer
nulls by permuting finite `final_spin_median` values across event labels while
keeping event identity, non-final-spin variables, missingness, and the historical
weighting/selector architecture fixed. Each outer catalog is passed through the
same 168 scan definitions.

Two independent-size calibrations were run:

| quantity | 200 catalogs | **1000 catalogs** |
|---|---:|---:|
| fixed-membership scans | 147 | 147 |
| dynamic-selector scans | 21 | 21 |
| MC floor | 0.004975 | **0.000999** |
| histogram chi-square global p | 0.094527 | **0.103896** |
| histogram chi-square global Z | 1.31338 sigma | **1.25966 sigma** |
| count-below-0.10 global p | 0.054726 | **0.062937** |

The 1000-catalog run used seed `20260901`. Its rank-calibrated observed p-vector
has an analytic independent-uniform conversion of `11.4562 sigma`, but that is
not a valid global significance because the p-values are discrete and strongly
dependent. The complete-vector empirical result is `p = 0.1038961`, or
`1.25966 sigma`, for this specific aggregate statistic and null.

### Prespecified null-sensitivity study

The unrestricted N0 result raised a legitimate question: does freely permuting
final spin across all masses destroy physically meaningful broad structure and
therefore over-correct the nominal 8.05-sigma base result?

To test that question without changing the observed statistic, the repo froze a
null-sensitivity protocol in [`GWTC4_NULL_SENSITIVITY.md`](GWTC4_NULL_SENSITIVITY.md)
before running three more models. The historical 168 scans, selectors, k-grid,
baseline degree, and observed scan-max statistics remained unchanged.

All sensitivity runs used **1000 complete outer catalogs**, and all first
reproduced the historical scan implementation with the same floating-point
errors above. The test suite passed **15/15** before execution.

| null | structure preserved | observed rank-p histogram `chi^2` | nominal analytic Z* | **empirical global p** | **empirical global Z** | count `p<0.10` global p |
|---|---|---:|---:|---:|---:|---:|
| N0 full event-label | non-spin catalog quantities + selector dependence | 167.119 | 11.456 sigma | **0.103896** | **1.260 sigma** | 0.062937 |
| N1 mass-stratified | broad total-mass association | 79.857 | 7.121 sigma | **0.344655** | **0.400 sigma** | 0.171828 |
| N2 mass+precision stratified | broad total-mass + final-spin precision association | 67.417 | 6.310 sigma | **0.428571** | **0.180 sigma** | 0.128871 |
| N3 smooth mass-spin residual | smooth degree-2 broad mass-spin relation | 181.583 | 12.040 sigma | **0.063936** | **1.523 sigma** | 0.037962 |

\*The nominal analytic Z values in this table are **diagnostic only** and are
not valid global significances for the dependent rank-calibrated vectors.

The primary empirical complete-vector result does not cross the prespecified
escalation threshold in any of the four null models. N3 is the most anomalous
primary result at `p = 0.063936` (`1.52255 sigma`); its secondary
count-below-0.10 statistic gives `p = 0.037962` (~`1.77 sigma`).

No structure-preserving null produced primary `p < 0.01`, and none approached
`p < 0.0027` (~3 sigma). Thus the nominal 8.0522-sigma conversion is not a
calibrated global significance for the historical p-vector statistic once the
catalog/selector dependence represented by N0–N3 is included.

That statement is deliberately narrower than saying the observed GWTC
structure has been "explained away." The nulls answer whether the **chosen
aggregate statistic** distinguishes the real catalog from catalogs preserving
specified correlation classes. A null that preserves broad mass-spin structure
may preserve structure of physical interest as well; causal attribution
requires a different test.

Full numerical outcomes and the prespecified interpretation are recorded in
[`GWTC4_NULL_SENSITIVITY_RESULTS.md`](GWTC4_NULL_SENSITIVITY_RESULTS.md).

### Historical GWTC-4 conclusion

The **8.0522-sigma historical result remains a real descriptive property of the
base scan ensemble under an independent-uniform reference**. It is useful and
should not be described as meaningless.

Across N0–N3, however, the same complete-vector histogram statistic has empirical
global significance of only approximately **0.18 sigma to 1.52 sigma**. The
proper conclusion is therefore:

> The historical 168-scan aggregate does **not currently establish** a globally
> calibrated >5-sigma GWTC-4 population anomaly under the tested dependence-aware
> nulls. Those null results do **not** establish a unique conventional physical
> explanation for the observed structure.

The newer V1/V3/V4 program addresses this limitation by freezing a specific
chirp-mass residual mode and asking whether its **frequency, amplitude, and
phase organization predict holdout data and survive non-periodic structured
population challenges**.

## Strict frozen holdout V1

The strict frozen holdout result is committed in
`tables/gwtc_frozen_holdout_result.csv`:

| quantity | result |
|---|---:|
| training N | 173 |
| holdout positive N | 104 |
| frozen `k` | `9.6023256203` |
| frozen amplitude | `0.4710666819` |
| holdout Delta deviance | `10.0354251` |
| null simulations | 10000 |
| holdout p | **0.0012998700** |
| verdict | **PASS_FIXED_MODE** |

This test is materially different from the historical p-vector aggregate: one
mode is selected on training data and then evaluated on a declared holdout
without rescanning the holdout.

## V3 — unbinned KDE fixed-mode robustness

V3 removes two major degrees of freedom from the earlier binned analysis:
histogram binning and polynomial background choice. The smooth baseline is a
Gaussian KDE selected by training-only leave-one-out likelihood, and the
residual is an unbinned normalized exponential tilt.

The committed V3 holdout result is:

| quantity | result |
|---|---:|
| selected KDE bandwidth | `0.1967599704` |
| frozen `k` | `9.7941176471` |
| frozen amplitude | `0.3850095492` |
| frozen phase | `1.3386452400 rad` |
| holdout Delta 2logL | **10.01220295** |
| null simulations | 10000 |
| null >= observed | **8** |
| holdout p | **0.0008999100** |
| verdict | **PASS_ROBUSTNESS_FIXED_MODE** |

Because GWTC-5 had already been inspected in V2, V3 is correctly classified as
a robustness result rather than an independent future-catalog replication.
Nevertheless, it establishes that the fixed-mode holdout effect is not solely a
consequence of the original histogram binning or degree-3 polynomial baseline.

See [V3_REPRODUCE.md](V3_REPRODUCE.md).

## V4 — frozen non-periodic structured-population challenge

V4 is the current attribution-focused test. It was frozen before its evaluator
is run. It fixes the externally published `k = 9.7` and training-fit residual
coefficients, including phase and amplitude, while the null population is a
**non-periodic** continuous broken power law plus truncated Gaussian peak. Three
predeclared selection-weight scenarios use `gamma = 0.0, 1.5, 2.5`.

The committed frozen V4 model has:

| quantity | frozen value |
|---|---:|
| `k` | `9.7` |
| amplitude | `0.3853480195` |
| phase | `1.0982350915 rad` |
| training Delta 2logL | `12.29726575` |
| equivalent mass ratio `exp(2pi/k)` | `1.9112377381` |

The V4 evaluator cannot scan `k`, refit `a,b`, refit the structured population,
or change the selection scenarios. Its result should be interpreted as a
**discrimination test between the exact frozen mode and the declared
non-periodic structured null family**, not as a binary judgment on whether
GWTC contains real structure.

At the time of this status update, no committed
`tables/gwtc_v4_structured_null_result.csv` is present. The V4 result is
therefore **pending**, not inferred from earlier null studies.

See [V4_REPRODUCE.md](V4_REPRODUCE.md).

## Catalog source and event counts

- **Source**: official GWOSC Event API,
  `https://gwosc.org/eventapi/json/allevents/` (see [PROVENANCE.md](PROVENANCE.md)).
- **Total event entries** across all catalogs: 671.
- **Unique events** (by common name): 433.
- **Cumulative deduplicated, `p_astro >= 0.5`**: **386 events** — consistent
  with the public GWTC-5.0 abstract's "roughly 390 transients".
- Catalogs include GWTC-4.0 (129), GWTC-4.1 (140), GWTC-5.0 (161), plus earlier
  releases.

## Original binned diagnostic and stress tests

The original current-harness primary used `M_chirp`, `ell = ln(M_chirp)`, a
simple polynomial Poisson baseline, and a scanned log-periodic residual. It
returned a small bootstrap p-value but was bin-fragile and coordinate
non-specific.

**Bin-count stress** (`tables/gwtc_bin_stress.csv`):

| B  | k_star | n_star | global_p |
|----|--------|--------|----------|
| 20 | 32.70  | 23.13  | 0.0025   |
| 30 |  4.48  |  3.17  | 0.0025   |
| 40 |  4.48  |  3.17  | 0.0025   |
| 50 |  4.48  |  3.17  | 0.0025   |

The coarse `B = 20` high-`k` solution is aliased; the finer bins settle near
`k ~ 4.48`, `n_star ~ 3.17`. The original scanned candidate therefore remains
**bin-fragile**.

**Variable stress** (`tables/gwtc_variable_stress.csv`): the bounded non-scale
coordinate `|chi_eff|` is also significant under the original scan, so that
scanned statistic is **not specific to scale-like variables**.

**Negative controls** (`tables/gwtc_control_results.csv`):

| control          | global_p range | reading |
|------------------|----------------|---------|
| `smooth_resample`| 0.19 – 0.61    | structureless → does not beat primary |
| `uniform_ell`    | 0.16 – 0.91    | structureless → does not beat primary |
| `jitter`         | 0.0025 – 0.0075| candidate survives posterior jitter |
| `bounded_chi_eff`| 0.0025         | non-scale coordinate is also significant |

### Interpretation of the original scanned candidate

The original binned scan demonstrates that the GWTC catalog contains structure
that a deliberately simple baseline does not absorb. Its bin fragility and
coordinate non-specificity mean that **this particular scanned statistic is not
sufficient to identify a universal WCT winding**.

The data do not, by themselves, establish that the residual is caused by the
known mass-function peaks, nor do they establish a WCT origin. Both conventional
population structure and additional phase-organized residual structure remain
hypotheses to discriminate. V3 and V4 were built specifically to move beyond
that ambiguity.

The historical aggregate, V1 fixed holdout, V3 unbinned robustness result, and
V4 structured-null challenge should therefore be reported separately rather
than compressed into one verdict label.
