# RESULTS

WCT-motivated GWTC log-domain residual diagnostic. **This is a diagnostic
harness only — it does not prove WCT and does not replace LVK population
inference.** All results below are reproducible via [REPRODUCE.md](REPRODUCE.md).

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
`1.25966 sigma`.

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

The primary empirical complete-vector result is nonsignificant in **all four
null models**. N1 and N2 are substantially *less* anomalous than N0. N3 is the
most favorable structure-preserving null for an anomaly, but still gives only
`p = 0.063936` (`1.52255 sigma`) on the primary population statistic. Its
secondary count-below-0.10 statistic gives `p = 0.037962` (~`1.77 sigma`), far
from the prespecified ~3-sigma escalation threshold and not an independent
replication.

No structure-preserving null produced primary `p < 0.01`, and none approached
`p < 0.0027` (~3 sigma). Thus the reduction of the nominal 8.0522-sigma base
aggregate is **not an artifact unique to unrestricted final-spin exchangeability**.
It persists when broad mass association, broad mass plus measurement precision,
and a smooth mass-spin relation are deliberately preserved.

Full numerical outcomes and the prespecified interpretation are recorded in
[`GWTC4_NULL_SENSITIVITY_RESULTS.md`](GWTC4_NULL_SENSITIVITY_RESULTS.md).

### Historical GWTC-4 conclusion

The **8.0522-sigma historical result remains a real descriptive property of the
base scan ensemble under an independent-uniform reference**. It is useful and
should not be described as meaningless. However, it does not survive as a
catalog-level physical anomaly when the complete dependent scan architecture is
calibrated.

Across the four tested dependence-aware nulls, the primary global significance
ranges from approximately **0.18 sigma to 1.52 sigma**. Therefore this repository
contains **no globally calibrated >5-sigma GWTC-4 population result** from the
historical 168-scan aggregate.

This conclusion does not prove that GWTC contains no physical structure. The
next materially different test would be a generative astrophysical population
null incorporating source-population structure, selection effects, and
parameter-estimation uncertainty. Such a model should be motivated by physical
fidelity rather than by whether it increases sigma.

## Catalog source and event counts

- **Source**: official GWOSC Event API,
  `https://gwosc.org/eventapi/json/allevents/` (see [PROVENANCE.md](PROVENANCE.md)).
- **Total event entries** across all catalogs: 671.
- **Unique events** (by common name): 433.
- **Cumulative deduplicated, `p_astro >= 0.5`**: **386 events** — consistent
  with the public GWTC-5.0 abstract's "roughly 390 transients".
- Catalogs include GWTC-4.0 (129), GWTC-4.1 (140), GWTC-5.0 (161), plus earlier
  releases.

## Variables available

`M_chirp`, `M_total`, `mass_1_source`, `mass_2_source`, `M_final`, `D_L`,
`redshift`, `chi_eff` (bounded, diagnostic only), `p_astro`, `far`.
`E_rad` and peak luminosity are **not** in the GWOSC summary table → reported
unavailable.

## Declared primary

- **Variable**: `M_chirp` (source-frame chirp mass), coordinate `ell = ln(M_chirp)`.
- **Subset**: cumulative deduplicated, `p_astro >= 0.5` (277 events with a
  finite positive `M_chirp`).
- **k-grid**: 0.5 .. 40.0, 120 points. **Null**: parametric Poisson bootstrap,
  `N_null = 1000`, seed 12345.

## Primary result (bin_count = 20)

| quantity      | value      |
|---------------|------------|
| `k_star`      | 32.70      |
| `n_star`      | 23.13      |
| `DeltaD_star` | 80.41      |
| `local_p`     | 3.5e-18 (diagnostic only) |
| `global_p`    | **0.0040** |
| `N_null`      | 1000       |

## Stress tests

**Bin-count stress** (`tables/gwtc_bin_stress.csv`):

| B  | k_star | n_star | global_p |
|----|--------|--------|----------|
| 20 | 32.70  | 23.13  | 0.0025   |
| 30 |  4.48  |  3.17  | 0.0025   |
| 40 |  4.48  |  3.17  | 0.0025   |
| 50 |  4.48  |  3.17  | 0.0025   |

`n_star` is **not stable** across the full declared bin set (CV = 1.06). The
coarse `B = 20` produces an **aliased high-k peak** (`k ~ 33`, above the
`B = 20` Nyquist frequency). For finer binning (`B >= 30`) the peak settles at
`k ~ 4.48`, `n_star ~ 3.17`. The harness therefore flags the candidate as
**bin-fragile**.

**Variable stress** (`tables/gwtc_variable_stress.csv`): all scale variables
return `global_p ~ 0.0025-0.005` at high, bin-fragile `k`. Notably the bounded
non-scale coordinate `|chi_eff|` is **equally significant** (`global_p ~ 0.0025`)
→ the detected mode is **not specific to scale-like variables**.

**Threshold / observing-run stress** (`tables/gwtc_threshold_stress.csv`):
significant under `p_astro >= 0.5`, `p_astro >= 0.9`, `FAR < 1/yr`, BBH-only,
O3a, O4a, O4b; weaker in O3b (`global_p ~ 0.07`).

**Negative controls** (`tables/gwtc_control_results.csv`):

| control          | global_p range | reading |
|------------------|----------------|---------|
| `smooth_resample`| 0.19 – 0.61    | structureless → **does not beat primary** (null well-calibrated) |
| `uniform_ell`    | 0.16 – 0.91    | structureless → **does not beat primary** |
| `jitter`         | 0.0025 – 0.0075| candidate **survives** posterior jitter (robustness positive) |
| `bounded_chi_eff`| 0.0025         | non-scale coordinate is **also significant** (coordinate non-specific) |

## Verdict

**VERDICT: PARTIAL — Reliability Class III** (`tables/gwtc_verdict.csv`,
`outputs/summary/VERDICT.txt`).

The primary is statistically significant under the 1000-replicate Poisson
bootstrap, and structureless controls do **not** beat it (the scan does not
manufacture significance on structureless data). However it fails two
robustness layers:

1. **Bin-fragility**: `n_star` is not stable across bin counts (coarse `B = 20`
   gives an aliased high-k peak; `B >= 30` collapses to `n_star ~ 3.17`).
2. **Coordinate non-specificity**: the bounded non-scale coordinate `|chi_eff|`
   is equally significant, so the mode is not specific to scale-like variables.

### Interpretation (bounded)

The large deviance improvement almost certainly reflects **known astrophysical
GWTC mass-function structure** (the well-documented ~10 M☉ and ~35 M☉ chirp/
component-mass over-densities) that a deliberately simple degree-3 Poisson
baseline cannot absorb — **not** a universal WCT log-periodic winding. This is
a **candidate log-domain residual** that is **bin-fragile and
coordinate-non-specific**, and it **requires independent confirmation** with a
more flexible baseline before any stronger statement. The harness remains
**diagnostic only**.

No result reaches **Class I**. The current best class is **Class III**. There
are no Class IV (pure-null) primary findings, but the structureless controls
behave as expected Class IV references.

![primary scan](figures/scan_M_chirp_B30.png)
