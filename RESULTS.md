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

This closes the earlier reconstruction gap: the historical selector/scanner
architecture is now reproducible from the restored files.

### Correlation-aware catalog-level outer null

`scripts/generate_gwtc4_population_null_matrix.py` generates whole-catalog outer
nulls by permuting finite `final_spin_median` values across event labels while
keeping event identity, non-final-spin variables, missingness, and the historical
weighting/selector architecture fixed. Each outer catalog is passed through the
same 168 scan definitions.

For the first **200 outer catalogs**:

- fixed-membership scans: 147;
- dynamic-selector scans rebuilt per permutation: 21;
- per-scan Monte Carlo floor: `1/201 = 0.0049751244`;
- observed outer-calibrated count `p < 0.05`: 43;
- observed outer-calibrated count `p < 0.10`: 52.

The outer-calibrated observed p-vector itself has
`chi^2 = 164.61905` and an analytic independent-uniform conversion of
`p = 3.6006e-30` (`11.3526 sigma`). **That analytic conversion is not a valid
global significance**, because these rank-calibrated p-values are both discrete
and strongly dependent.

The valid comparison is against the complete correlated null-catalog
population statistics:

| end-to-end statistic | empirical result |
|---|---:|
| histogram chi-square global p | **0.094527363** |
| histogram chi-square global Z | **1.31338 sigma** |
| count-below-0.10 global p | **0.054726368** |
| outer catalogs | 200 |
| MC resolution floor | 0.0049751244 |

**Result:** under the stated final-spin event-label permutation null, the
historical GWTC-4 population anomaly is **not globally significant**. The
nominal 8.0522-sigma historical value collapses to an empirical
**1.31-sigma catalog-level result** after selector/event-overlap dependence is
included.

Therefore this repository contains **no globally calibrated >5-sigma GWTC-4
population result**. The old >5-sigma values are retained only to document the
historical marginal calculation and the size of the dependence correction.

This permutation null is not a complete LVK astrophysical population model; it
tests a specific exchangeability null for final spin under the recovered
selector architecture. Additional astrophysical/measurement-aware nulls are
useful robustness checks, but the current correlated result already invalidates
the earlier independent-uniform >5-sigma interpretation.

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
