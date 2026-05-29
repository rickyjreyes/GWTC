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
