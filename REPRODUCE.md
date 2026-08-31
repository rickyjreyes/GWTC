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

## 10. Historical GWTC-4 168-scan population calibration

The historical analysis reported 168 scan-max p-values. The analytic population
statistic is reproduced from the observed p-value vector with:

```bash
python scripts/run_gwtc_population_global_null.py \
  --observed tables/gwtc4_population_observed.csv \
  --output tables/gwtc4_population_global_result.json
```

For the **actual correlated catalog-level global-null test**, supply a matrix
created by rerunning the identical 168-scan workflow on complete null catalogs:

```bash
python scripts/run_gwtc_population_global_null.py \
  --observed tables/gwtc4_population_observed.csv \
  --null-matrix tables/gwtc4_population_null_matrix.csv \
  --output tables/gwtc4_population_global_result.json
```

Required formats:

- `gwtc4_population_observed.csv`: exactly 168 rows with one p-value column
  named `p_scanmax`, `global_p`, `p_value`, or `p` (or specify `--p-column`).
- `gwtc4_population_null_matrix.csv`: one complete null catalog per row and
  exactly 168 p-value columns in the same scan-definition order as the observed
  vector. Metadata columns may be present if they are nonnumeric; alternatively
  prefix all scan columns (for example `scan_000` ... `scan_167`) and pass
  `--null-prefix scan_`.

The historical statistic uses bin edges
`[0,.05,.10,.20,.30,.40,.50,.60,.70,.80,.90,1.0]`. This is the 11-bin
construction consistent with the reported expected counts 8.4, 8.4, then 16.8
for `N=168`.

**Do not use 168 independent Uniform(0,1) draws per row as the null matrix.**
That would erase the overlap/correlation among ranking and subset selectors and
would only reproduce the nominal analytic test. To support a globally
calibrated discovery-level statement, each null row must come from one complete
synthetic/null catalog passed through the same 168-scan pipeline.

The empirical Monte Carlo p-value uses the standard plus-one correction
`(1 + N_ge) / (1 + N_null_catalogs)`. Therefore direct resolution below the
one-sided 5-sigma threshold requires roughly 3.49 million null catalogs if no
null realization exceeds the observed aggregate statistic.

## Outputs

- `tables/gwtc_*.csv` — canonical results (one schema; see METHOD.md §6).
- `outputs/summary/VERDICT.txt` — human-readable verdict.
- `outputs/summary/scan_*.csv`, `outputs/summary/*.png` — scan curves / plots.
- `tables/gwtc4_population_global_result.json` — analytic and, when supplied,
  empirical catalog-level population calibration.

## Seeds

Default seeds are explicit: primary/stress `--seed 12345`, controls
`--seed 2024`. Change them to test stability; the verdict should not hinge on a
single seed.
