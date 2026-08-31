# REPRODUCE

Every step is deterministic given the explicit random seeds. Every script
supports `--help`.

## 0. Environment

```bash
pip install -r requirements.txt
```

Requires: numpy, pandas, scipy, matplotlib, pytest. No seaborn.

## 1. Fetch official catalog + provenance

```bash
python scripts/fetch_gwtc_catalog.py
# -> data/gwtc_allevents_raw.json, data/manifest.csv (URL, access date, sha256)
```

## 2. Inspect schema / variable availability

```bash
python scripts/inspect_catalog_schema.py
```

## 3. Build the clean per-event table

```bash
python scripts/build_gwtc_table.py
# -> tables/gwtc_events_clean.csv
```

## 4. Smoke run (fast)

```bash
python scripts/run_gwtc_log_scan.py --variable M_chirp --cumulative \
    --p-astro-min 0.5 --bins 20 --null-n 100 --k-min 0.5 --k-max 40 --n-k 60
pytest -q
```

## 5. Full primary run + null (PASS-grade null_n >= 1000)

```bash
python scripts/run_gwtc_log_scan.py --variable M_chirp --cumulative \
    --p-astro-min 0.5 --bins 20 --null-n 1000 --n-k 120 --seed 12345
python scripts/run_gwtc_nulls.py --variable M_chirp --cumulative \
    --p-astro-min 0.5 --bins 20 --null-n 1000 --n-k 120 --seed 12345
```

## 6. Stress tests

```bash
# bin-count stress: 20, 30, 40, 50
python scripts/run_gwtc_bin_stress.py --variable M_chirp --cumulative \
    --p-astro-min 0.5 --null-n 400 --n-k 120

# variable stress: ln(M_chirp), ln(M_total), ln(D_L), ln(M_final), ln(redshift),
# plus bounded diagnostic |chi_eff|
python scripts/run_gwtc_variable_stress.py --cumulative --p-astro-min 0.5 \
    --bins 20 --null-n 400 --n-k 120 --include-bounded

# threshold + observing-run stress
python scripts/run_gwtc_threshold_stress.py --variable M_chirp \
    --bins 20 --null-n 400 --n-k 120
```

## 7. Negative controls

```bash
python scripts/run_gwtc_controls.py --variable M_chirp --cumulative \
    --p-astro-min 0.5 --bins 20 --null-n 400 --n-k 120
```

## 8. Master table + verdict

```bash
python scripts/make_gwtc_master_table.py
python scripts/make_gwtc_verdict.py --primary-variable M_chirp \
    --primary-subset cumulative_pa0.5
# -> tables/gwtc_verdict.csv, outputs/summary/VERDICT.txt
```

## 9. One-shot current diagnostic

```bash
bash run_all.sh    # convenience wrapper for steps 5-8 (see repo root)
```

## 10. Historical GWTC-4 168-scan population audit and calibration

### 10.1 Audit the reported aggregate arithmetic

This command requires no historical p-value file. It checks the reported
`N=168`, `chi^2=95.6`, 16 scans below `0.05`, and 47 scans below `0.10`:

```bash
python scripts/run_gwtc_population_global_null.py --audit-reported
```

The current audit should report approximately:

```text
chi2 = 95.6, df = 10
analytic chi2 p = 4.12e-16
nominal one-sided Z = 8.05 sigma

47 / 168 below 0.10:
  Poisson(lambda=16.8) upper tail = 1.16e-9   (~5.97 sigma)
  exact Binomial(168,0.10) upper tail = 4.85e-11 (~6.47 sigma)

16 / 168 below 0.05:
  Poisson(lambda=8.4) upper tail = 0.0125
  exact Binomial(168,0.05) upper tail = 0.0105
```

This exposes two historical labeling/arithmetic issues: the reported
`47/168` value labeled as a Poisson tail `<1e-10` is not the literal Poisson
result, although the exact binomial tail is `<1e-10`; and the reported
`16/168` Poisson-tail `1.4e-3` is not reproduced by the stated Poisson or exact
binomial calculation.

The paper text also says "11 equal bins," but its figure and expected counts
correspond to edges
`[0,.05,.10,.20,.30,.40,.50,.60,.70,.80,.90,1.0]`. The audit uses these
figure-compatible edges.

### 10.2 Recompute from the actual observed 168 p-values

If the historical observed vector is restored or regenerated:

```bash
python scripts/run_gwtc_population_global_null.py \
  --observed tables/gwtc4_population_observed.csv \
  --output tables/gwtc4_population_global_result.json
```

Required observed format:

- exactly 168 rows;
- one p-value column named `p_scanmax`, `global_p`, `p_value`, or `p`, or select
  the column explicitly with `--p-column`.

The exact historical 168-value vector and its selector manifest are **not
currently committed**, so this stage cannot yet reconstruct the published
histogram from repository data alone.

### 10.3 Actual correlated catalog-level global-null test

For the global test that matters scientifically, supply a matrix created by
rerunning the **identical 168-scan workflow on complete null catalogs**:

```bash
python scripts/run_gwtc_population_global_null.py \
  --observed tables/gwtc4_population_observed.csv \
  --null-matrix tables/gwtc4_population_null_matrix.csv \
  --output tables/gwtc4_population_global_result.json
```

Required null format:

- one complete null catalog per row;
- exactly 168 p-value columns in the same scan-definition order as the observed
  vector;
- metadata columns may be present if they are nonnumeric;
- alternatively prefix scan columns (for example `scan_000` ... `scan_167`)
  and pass `--null-prefix scan_`.

**Do not use 168 independent Uniform(0,1) draws per row as the null matrix.**
That erases the overlap/correlation among ranking and subset selectors and only
reproduces an independence-based analytic reference, not the actual global
null of the analysis workflow.

The empirical Monte Carlo p-value uses the plus-one correction
`(1 + N_ge) / (1 + N_null_catalogs)`. Directly resolving a one-sided 5-sigma
p-value with zero exceedances therefore requires roughly **3.49 million** null
catalogs. A smaller ensemble is still useful for testing whether the nominal
8-sigma-scale statistic collapses once scan dependence is included.

## Outputs

- `tables/gwtc_*.csv` — canonical results (one schema; see METHOD.md §6).
- `outputs/summary/VERDICT.txt` — human-readable verdict.
- `outputs/summary/scan_*.csv`, `outputs/summary/*.png` — scan curves / plots.
- `tables/gwtc4_population_global_result.json` — arithmetic audit and, when
  supplied, observed-vector and empirical catalog-level population calibration.

## Seeds

Default seeds are explicit: primary/stress `--seed 12345`, controls
`--seed 2024`. Change them to test stability; the verdict should not hinge on a
single seed.
