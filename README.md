# WCT-Motivated GWTC Log-Domain Residual Diagnostic

This repository tests **WCT-style log-domain residual structure** in
LIGO-Virgo-KAGRA Gravitational-Wave Transient Catalog (GWTC) event-catalog
variables.

> **Scope and honesty statement**
>
> - It **does not prove WCT**.
> - It **does not replace LVK population-inference analyses**.
> - It **does not claim a new gravitational-wave detection**.
> - It tests whether selected catalog-scale variables contain reproducible,
>   predictive, or otherwise nontrivial log-domain residual structure under
>   increasingly strict null models and holdout protocols.
>
> A null result does not by itself prove that the underlying catalog is
> structureless, and a positive result does not by itself identify the physical
> mechanism. The repository separates **descriptive structure**, **predictive
> holdout evidence**, and **source-attribution / null discrimination**.

## Current evidence hierarchy

The repository now contains several distinct GWTC tests that should not be
collapsed into a single sigma value.

| layer | question | current result | interpretation |
|---|---|---|---|
| Historical GWTC-4 aggregate | Are the 168 historical scan outputs unusually nonuniform under an independent-uniform reference? | nominal `8.0522 sigma` | **Yes descriptively**, but the 168 scans are strongly dependent. |
| GWTC-4 dependence-aware outer nulls | Is that same 168-scan p-vector rare after preserving selector/catalog dependence? | empirical global `0.18–1.52 sigma` across N0–N3 | The historical aggregate is **not globally calibrated >5 sigma** under these nulls. This does not prove the structure is nonphysical. |
| Strict frozen holdout V1 | Does a training-frozen chirp-mass mode predict the declared holdout? | `p = 0.00129987`, `PASS_FIXED_MODE` | **Positive fixed-mode holdout evidence.** |
| Unbinned KDE V3 | Does the frozen mode survive removal of histogram binning and replacement of the polynomial baseline by a training-only Gaussian KDE? | `Delta 2logL = 10.0122`, `8/10000` exceedances, `p = 0.00089991`, `PASS_ROBUSTNESS_FIXED_MODE` | **Positive unbinned fixed-mode robustness evidence.** GWTC-5 had already been inspected in V2, so this is robustness rather than a new independent catalog replication. |
| Structured astrophysical-null V4 | Can a training-fitted **non-periodic** broken-power-law + Gaussian-peak population with predeclared selection weighting reproduce the exact frozen `k = 9.7`, amplitude, and phase statistic? | **FROZEN; evaluation pending** | This is the current source-discrimination challenge. The null deliberately preserves broad population structure without inserting a periodic term. |

The strict holdout result is committed in
`tables/gwtc_frozen_holdout_result.csv`. The V3 result is committed in
`tables/gwtc_v3_unbinned_kde_holdout_result.csv`. The V4 frozen definition is
committed in `tables/gwtc_v4_frozen_structured_null.json`; see
[V4_REPRODUCE.md](V4_REPRODUCE.md) for the one-shot evaluation protocol.

The key interpretive distinction is:

> **Failure to distinguish an observed statistic from a structure-preserving
> null is not equivalent to proving that the null's physical mechanism caused
> the observation.** It means that statistic is not diagnostic between the
> observed catalog and that declared null family.

Conversely, a V4 pass would show that the exact frozen phase/amplitude/frequency
mode is uncommon under the declared non-periodic structured population family;
it would not by itself establish WCT or fully model LVK selection effects.

## GWTC-4 population-level aggregate result

The restored historical GWTC-4 subset analysis contains exactly **168 scans =
56 subset selectors x 3 final-spin pairs**. Recomputing the historical aggregate
from the original 168-row scan table gives a very large **nominal** departure
from an independent-uniform reference:

| statistic | result | interpretation |
|---|---:|---|
| historical p-vector histogram | `chi^2 = 95.630952`, `p = 4.06578e-16` | nominal `Z = 8.0522 sigma` |
| historical count `p < 0.10` | 47 / 168 | 5.97 sigma Poisson / 6.47 sigma exact-binomial marginally |
| historical count `p < 0.05` | 16 / 168 | about 2.24-2.31 sigma marginally |

The historical paper text says "11 equal bins," while its plotted bins and
reported expected counts correspond to `[0,.05]`, `[.05,.10]`, followed by
nine width-`.10` bins. The audit uses those figure-compatible bin edges because
they reproduce the stated expectations 8.4, 8.4, then 16.8 for `N = 168`.

### Reproduction and correlation-aware calibration

The exact historical scan implementation was restored and independently
reproduced before null generation:

- 86 events;
- 168 scans;
- 56 unique subsets;
- 3 final-spin scan pairs;
- maximum reconstructed `k_best` difference `7.11e-15`;
- maximum reconstructed `Delta chi^2` difference `7.43e-09`;
- **observed reproduction: PASS**.

The first full event-label permutation null (N0) was run at both 200 and 1000
catalogs. The independent 1000-catalog run is the primary N0 calibration:

| N0 statistic | 200 catalogs | **1000 catalogs** |
|---|---:|---:|
| population histogram chi-square global p | 0.094527 | **0.103896** |
| population histogram global Z | 1.313 sigma | **1.260 sigma** |
| excess count below `p < 0.10` global p | 0.054726 | **0.062937** |
| MC resolution floor | 0.004975 | **0.000999** |

The 1000-catalog rank-calibrated observed p-vector itself gives a nominal
independent-uniform `Z = 11.4562 sigma`, but that conversion is **not a valid
global significance** because those p-values are discrete and strongly
dependent. The complete-vector empirical comparison is the scientifically
relevant calibration for that specific aggregate statistic.

### Prespecified null-sensitivity study

To test whether unrestricted final-spin exchangeability was itself too
destructive, the observed 168-scan analysis was frozen and three additional
nulls were declared **before seeing their results**. Only the catalog-level null
model changed.

| null | structure deliberately preserved | empirical histogram global p | global Z | count `p<0.10` global p |
|---|---|---:|---:|---:|
| N0 full event-label permutation | selector/event overlap and non-spin catalog quantities | `0.103896` | `1.260 sigma` | `0.062937` |
| N1 mass-stratified permutation | broad total-mass association | `0.344655` | `0.400 sigma` | `0.171828` |
| N2 mass+precision stratified | broad mass and spin-precision association | `0.428571` | `0.180 sigma` | `0.128871` |
| N3 smooth mass-spin residual permutation | smooth degree-2 broad mass-spin relation | `0.063936` | `1.523 sigma` | `0.037962` |

All alternative runs used 1000 complete outer catalogs and reproduced the
historical scanner before scoring. None of the three structure-preserving
nulls produced `p < 0.01` on the primary complete-vector histogram statistic.
N3 is the most anomalous primary result at `1.52 sigma`; its secondary count
statistic is `p = 0.03796` (~`1.77 sigma`).

These results establish a limitation of the **historical aggregate statistic**:
its nominal 8.0522-sigma independent-uniform conversion is not a calibrated
catalog-level significance once the dependence architecture is included. They
do **not** establish that the underlying catalog structure has been causally
explained by those permutation models. A null that deliberately preserves a
class of correlations can also preserve structure of physical interest; the
question is whether the chosen statistic distinguishes the observed catalog
from that null family.

**Current historical conclusion:** the `8.0522 sigma` base-run value remains a
real descriptive statement about the historical scan-output distribution under
an independent-uniform reference. It is not meaningless and is retained. The
same p-vector histogram statistic, however, does **not currently establish a
globally calibrated >5-sigma GWTC-4 population anomaly** under N0–N3.

This is why the newer V1/V3/V4 line freezes a specific chirp-mass residual mode
and tests prediction and source discrimination directly rather than relying on
the aggregate distribution of 168 marginal p-values.

See [`GWTC4_NULL_SENSITIVITY.md`](GWTC4_NULL_SENSITIVITY.md) for the frozen
pre-result methodology and
[`GWTC4_NULL_SENSITIVITY_RESULTS.md`](GWTC4_NULL_SENSITIVITY_RESULTS.md) for the
recorded outcomes.

Reproduce the primary 1000-catalog N0 calibration:

```bash
python scripts/generate_gwtc4_population_null_matrix.py \
  --verify-only \
  --metadata tables/gwtc4_population_null_verify.json

python scripts/generate_gwtc4_population_null_matrix.py \
  --outer-n 1000 \
  --seed 20260901 \
  --metadata tables/gwtc4_population_null_metadata_1000.json \
  --output-matrix tables/gwtc4_population_null_matrix_1000.csv \
  --output-observed tables/gwtc4_population_observed_outercal_1000.csv \
  --output-delta tables/gwtc4_population_null_delta_matrix_1000.csv

python scripts/run_gwtc_population_global_null.py \
  --observed tables/gwtc4_population_observed_outercal_1000.csv \
  --p-column outer_global_p \
  --null-matrix tables/gwtc4_population_null_matrix_1000.csv \
  --output tables/gwtc4_population_global_result_1000.json
```

The exact historical scanner is preserved as
`scripts/gwtc4_wct_subset_scan_compound.py`, and the restored scan vector is
committed as `tables/gwtc4_population_observed.csv`. See
[REPRODUCE.md](REPRODUCE.md) for the full procedure and null-model limitations.

## Core hypothesis

If WCT-style curvature/winding structure appears in gravitational-wave
catalogs, it should appear as **stable residual phase organization in physically
scale-like logarithmic variables**, not merely as arbitrary small p-values.

For a positive, scale-like catalog variable `z` (e.g. source-frame chirp mass),
define the log coordinate `ell = ln(z)`. The basic residual family is

```text
A cos(k ell - phi)
```

or equivalently

```text
a cos(k ell) + b sin(k ell).
```

Exploratory scans may search `k`, but the stronger tests freeze `k`, amplitude,
phase, baseline definition, and data split before evaluation. The active-domain
winding number is `n_star = k_star * Delta_ell_A / (2*pi)`.

The current evidence hierarchy therefore prioritizes, in order:

1. independent/frozen prediction;
2. survival under non-polynomial and unbinned baselines;
3. survival under structured **non-periodic** population nulls;
4. only then conversion to stronger physical claims.

See [METHOD.md](METHOD.md), [V3_REPRODUCE.md](V3_REPRODUCE.md), and
[V4_REPRODUCE.md](V4_REPRODUCE.md) for the current protocols.

## Data source

Official **GWOSC Event API** (LIGO-Virgo-KAGRA / GWOSC):
`https://gwosc.org/eventapi/json/allevents/`. We use the official catalog
metadata rather than scraped secondary summaries. Provenance (URL, access
date, SHA-256, event counts) is recorded in [PROVENANCE.md](PROVENANCE.md)
and `data/manifest.csv`. Catalogs include GWTC-4.0, GWTC-4.1 and GWTC-5.0; the
cumulative deduplicated set with `p_astro >= 0.5` is ~386 events, consistent
with the public GWTC-5.0 abstract's "roughly 390 transients".

## Admissible variables

Only physically scale-like variables are used as **primary** evidence:

- `M_chirp`  — source-frame chirp mass (primary recommended)
- `M_total`  — source-frame total mass
- `D_L`      — luminosity distance
- `M_final`  — source-frame final mass (where available)
- `redshift` — where available

Bounded ratio/spin variables (e.g. `|chi_eff|`) are treated as **bounded,
non-scale, diagnostic-only** coordinates and can never produce a primary PASS.

## Quick start

```bash
pip install -r requirements.txt

# 1. Fetch official catalog metadata + record provenance
python scripts/fetch_gwtc_catalog.py

# 2. Inspect schema / variable availability
python scripts/inspect_catalog_schema.py

# 3. Build the clean per-event table
python scripts/build_gwtc_table.py

# 4. Run tests
pytest -q
```

For the historical diagnostic sequence see [REPRODUCE.md](REPRODUCE.md).
For the stronger frozen fixed-mode sequences see
[V2_REPRODUCE.md](V2_REPRODUCE.md), [V3_REPRODUCE.md](V3_REPRODUCE.md), and
[V4_REPRODUCE.md](V4_REPRODUCE.md).

## Repository layout

```text
scripts/   historical scanners, frozen holdout evaluators, and structured nulls
tables/    canonical results, frozen model artifacts, and manifests
data/      raw catalog payload + provenance manifest
outputs/   diagnostic scan curves and summaries
tests/     pipeline correctness, freeze guards, null tests, overclaim guards
```
