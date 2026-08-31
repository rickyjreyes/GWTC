# WCT-Motivated GWTC Log-Domain Residual Diagnostic

This repository tests **WCT-style log-domain residual structure** in LIGO-Virgo-KAGRA Gravitational-Wave Transient Catalog (GWTC) event-catalog variables.

> **Scope and honesty statement**
>
> - It **does not prove WCT**.
> - It **does not replace LVK population-inference analyses**.
> - It **does not claim a new gravitational-wave detection**.
> - Positive branches are treated as a **candidate log-domain residual** requiring independent confirmation and stronger source discrimination.
>
> The repository separates **descriptive structure**, **predictive holdout evidence**, and **source/null discrimination**. A null-model result is evidence about a statistic under that declared null family; it is not automatically a causal explanation of the observed catalog.

## Current evidence hierarchy

| layer | current result | interpretation |
|---|---|---|
| Historical GWTC-4 aggregate | nominal `8.0522 sigma` under an independent-uniform p-vector reference | Real descriptive nonuniformity, but the 168 scans are strongly dependent. |
| GWTC-4 dependence-aware N0-N3 nulls | empirical global `0.18-1.52 sigma` for the historical aggregate statistic | The old 168-scan aggregate is not globally calibrated >5 sigma under those nulls; that does not establish a unique physical explanation. |
| Strict frozen holdout V1 | `p = 0.00129987`, `PASS_FIXED_MODE` | Positive fixed-mode holdout evidence. |
| Unbinned KDE V3 | `Delta 2logL = 10.0122`, `8/10000`, `p = 0.00089991`, `PASS_ROBUSTNESS_FIXED_MODE` | The frozen mode survives removal of histogram binning and replacement of the polynomial baseline by a training-only KDE. |
| Structured non-periodic population-null V4 | **`PASS_ALL_STRUCTURED_NULLS`**; worst-case `p = 0.00559944` | The exact frozen `k = 9.7`, amplitude and phase statistic is uncommon under all three declared non-periodic structured-population scenarios. |

## V4 result

V4 was frozen before evaluation. The evaluator could not scan `k`, refit the residual coefficients, refit the structured population, or change the predeclared selection scenarios on the holdout.

Frozen signal:

- variable: `M_chirp`
- coordinate: `ell = log(M_chirp)`
- `k = 9.7`
- amplitude `= 0.3853480195`
- phase `= 1.0982350915 rad`
- holdout N `= 104`
- observed holdout `Delta 2logL = 9.95018558`

Declared non-periodic null family: continuous broken power law + truncated Gaussian peak, with selection weighting `(Mchirp/Mchirp_max)^gamma` for `gamma = 0.0, 1.5, 2.5`.

| gamma | null >= observed | empirical p | verdict |
|---:|---:|---:|---|
| `0.0` | `55/10000` | `0.0055994401` | `PASS_STRUCTURED_NULL_SCENARIO` |
| `1.5` | `1/10000` | `0.0001999800` | `PASS_STRUCTURED_NULL_SCENARIO` |
| `2.5` | `2/10000` | `0.0002999700` | `PASS_STRUCTURED_NULL_SCENARIO` |

Overall verdict: **`PASS_ALL_STRUCTURED_NULLS`**.

This is a robustness/source-discrimination result, not an independent future-catalog replication, because GWTC-5 had already been inspected in earlier V2/V3 work. The V4 null is also phenomenological rather than a full LVK hierarchical population model with event-posterior propagation and injection-calibrated selection effects.

See:

- `GWTC_V4_STRUCTURED_NULL_RESULT.md`
- `tables/gwtc_v4_frozen_structured_null.json`
- `tables/gwtc_v4_structured_null_result.csv`
- `V4_REPRODUCE.md`

## Why the historical 8.0522 sigma is kept separate

The restored historical analysis contains **168 scans = 56 subset selectors x 3 final-spin pairs**. Its p-vector is strongly nonuniform under an independent-uniform reference, giving the literal nominal `8.0522 sigma` calculation. That number is retained as a descriptive property of the historical scan ensemble.

However, N0-N3 whole-catalog dependence-aware calibrations give only about `0.18-1.52 sigma` for that particular aggregate statistic. The correct conclusion is therefore narrow: the nominal 8.0522-sigma conversion is not a calibrated catalog-level significance for the dependent 168-scan architecture. It does **not** follow that the underlying structure has been causally explained away.

The newer V1/V3/V4 line was built to ask a different question: can a **specific frozen phase/frequency/amplitude mode** predict holdout data and survive progressively stronger non-periodic nulls?

## Core hypothesis

If WCT-style curvature/winding structure appears in gravitational-wave catalogs, it should appear as **stable residual phase organization in physically scale-like variables**, not merely as arbitrary small p-values.

For a positive scale-like variable `z`, define `ell = ln(z)` and test a residual family

```text
A cos(k ell - phi)
```

or equivalently

```text
a cos(k ell) + b sin(k ell).
```

Exploratory scans may search `k`; stronger tests freeze `k`, amplitude, phase, baseline definition, and data split before evaluation.

## Controls and interpretation

The repo retains negative control results, failed robustness branches, historical diagnostics, and model limitations. Bounded/non-scale variables such as `|chi_eff|` are diagnostic-only and cannot produce a primary physical PASS.

The strongest missing conventional challenge is now a **posterior-aware hierarchical population/selection test** using detector injections or an equivalent validated selection model while keeping the residual mode externally frozen. The strongest future evidentiary test is an independent future GWTC release evaluated once against a prediction frozen before release.

## Data and reproduction

Source: official GWOSC Event API. Provenance is recorded in `PROVENANCE.md` and `data/manifest.csv`.

Key files:

- `RESULTS.md` — detailed evidence hierarchy and numerical results
- `METHOD.md` — definitions and verdict rules
- `REPRODUCE.md` — historical diagnostic reproduction
- `V3_REPRODUCE.md` — unbinned KDE fixed-mode robustness
- `V4_REPRODUCE.md` — frozen structured-population-null challenge
- `GWTC_V4_STRUCTURED_NULL_RESULT.md` — V4 result record
- `TODO.md` — remaining tests

Run the test suite with:

```bash
python -m pytest -q --basetemp=.pytest_tmp
```

This harness is deliberately designed to retain FAIL/PARTIAL outcomes and **requires independent confirmation** before stronger physical claims.