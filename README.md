# WCT-Motivated GWTC Log-Domain Residual Diagnostic

This repository tests **WCT-style log-domain residual structure** in
LIGO-Virgo-KAGRA Gravitational-Wave Transient Catalog (GWTC) event-catalog
variables.

> **Scope and honesty statement**
>
> - It **does not prove WCT**.
> - It **does not replace LVK population-inference analyses**.
> - It **does not claim a new gravitational-wave detection**.
> - It only tests whether selected catalog-scale variables contain stable
>   log-domain residual modes after null testing and stress checks.
>
> This is a **WCT-motivated diagnostic** only. A negative result does not
> disprove WCT, and a positive result is a **candidate log-domain residual**
> that **requires independent confirmation**. The harness is designed to be
> able to return **PASS, PARTIAL, FAIL, or INCOMPLETE**, and negative /
> negative-control results are retained as valuable.

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

### Correlation-aware outer-null result

The exact historical scan implementation was restored and independently
reproduced before null generation:

- 86 events;
- 168 scans;
- 56 unique subsets;
- 3 final-spin scan pairs;
- maximum reconstructed `k_best` difference `7.11e-15`;
- maximum reconstructed `Delta chi^2` difference `7.43e-09`;
- **observed reproduction: PASS**.

A first 200-catalog outer permutation ensemble gave a population histogram
`p = 0.094527` (`Z = 1.313 sigma`). A separate **1000-catalog** run with a new
seed (`20260901`) reproduced the same conclusion and is the primary quoted
correlation-aware calibration:

| correlation-aware statistic | 200 catalogs | **1000 catalogs** |
|---|---:|---:|
| population histogram chi-square global p | 0.094527 | **0.103896** |
| population histogram global Z | 1.313 sigma | **1.260 sigma** |
| excess count below `p < 0.10` global p | 0.054726 | **0.062937** |
| observed count below `p < 0.05` | 43 | **42** |
| observed count below `p < 0.10` | 52 | **52** |
| MC resolution floor | 0.004975 | **0.000999** |

For the 1000-catalog run, the outer-calibrated observed p-vector itself has
`chi^2 = 167.11905` and an analytic independent-uniform conversion of
`p = 1.09489e-30` (`11.4562 sigma`). **That analytic conversion is not a valid
global significance**, because these rank-calibrated p-values are discrete and
strongly dependent. The scientifically relevant result is the empirical
complete-vector comparison above: **`p = 0.103896`, `Z = 1.25966 sigma`**.

**Conclusion:** under this stated final-spin event-label permutation null, the
historical population anomaly is **not globally significant**. The earlier
8.0522-sigma historical value is retained only as a nominal independence-based
diagnostic; it must not be described as a >5-sigma catalog-level result. The
200- and 1000-catalog calibrations agree closely, showing that the dependence
among overlapping subsets/selectors explains the apparent population excess
under this null.

This outer null is not a full LVK astrophysical population model, so the result
is a statement about this explicit permutation null rather than a universal
proof of no structure. Additional null models can test robustness, but the
current repository contains **no globally calibrated >5-sigma GWTC-4 population
result**.

Reproduce the primary 1000-catalog robustness calibration:

```bash
# Verify restored historical scan implementation
python scripts/generate_gwtc4_population_null_matrix.py \
  --verify-only \
  --metadata tables/gwtc4_population_null_verify.json

# Generate the 1000-catalog correlated outer null
python scripts/generate_gwtc4_population_null_matrix.py \
  --outer-n 1000 \
  --seed 20260901 \
  --metadata tables/gwtc4_population_null_metadata_1000.json \
  --output-matrix tables/gwtc4_population_null_matrix_1000.csv \
  --output-observed tables/gwtc4_population_observed_outercal_1000.csv \
  --output-delta tables/gwtc4_population_null_delta_matrix_1000.csv

# Calibrate the complete population statistic
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

This historical aggregate result is intentionally kept separate from the
current stricter diagnostic verdict in [RESULTS.md](RESULTS.md). The current
harness verdict remains **PARTIAL — Reliability Class III**.

## Core hypothesis

If WCT-style curvature/winding structure appears in gravitational-wave
catalogs, it should appear as **stable residual structure in physically
scale-like logarithmic variables**, not arbitrary transforms.

For a positive, scale-like catalog variable `z` (e.g. source-frame chirp mass),
define the log coordinate `ell = ln(z)`, bin it into an event-density field,
fit a smooth Poisson baseline `mu0(ell)`, and test a single log-periodic
residual mode:

```
log mu(ell; k, a, b, c) = log mu0(ell) + c + a cos(k ell) + b sin(k ell)
```

scanning `k` over a declared grid. The active-domain winding number is
`n_star = k_star * Delta_ell_A / (2*pi)`. Significance is judged by a
**parametric Poisson-bootstrap global p-value**, never by a local chi^2
p-value alone.

See [METHOD.md](METHOD.md) for full definitions and the PASS/PARTIAL/FAIL/
INCOMPLETE rules, and [RESULTS.md](RESULTS.md) for the current verdict.

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

# 4. Smoke run on ln(M_chirp)
python scripts/run_gwtc_log_scan.py --variable M_chirp --cumulative \
    --p-astro-min 0.5 --bins 20 --null-n 100 --k-min 0.5 --k-max 40

# 5. Tests
pytest -q
```

For the full reproduction sequence (all stress tests, controls, verdict) see
[REPRODUCE.md](REPRODUCE.md).

## Repository layout

```
scripts/   pipeline (fetch, build, scan, nulls, stress, controls, verdict)
tables/    canonical CSV outputs (one schema, see METHOD.md)
data/      raw catalog payload + provenance manifest
outputs/   per-k scan curves and the human-readable VERDICT.txt
tests/     pipeline correctness + overclaim guard
```

Every result CSV carries the same schema columns: `catalog_version, variable,
subset, bin_count, k_best, n_star, DeltaD_star, local_p, global_p, null_n,
verdict_label, notes`.
