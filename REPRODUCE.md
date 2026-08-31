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

The recovered historical run contains exactly **168 scans = 56 subset
selectors x 3 final-spin pairs**. The original 168-row scan result table is
committed as `tables/gwtc4_population_observed.csv`, and the historical scanner
is preserved as `scripts/gwtc4_wct_subset_scan_compound.py`.

### 10.1 Audit the reported aggregate arithmetic

```bash
python scripts/run_gwtc_population_global_null.py --audit-reported
```

The arithmetic audit reports approximately:

```text
chi2 = 95.6, df = 10
analytic chi2 p = 4.12e-16
nominal one-sided Z = 8.05 sigma

47 / 168 below 0.10:
  Poisson(lambda=16.8) upper tail = 1.16e-9      (~5.97 sigma)
  exact Binomial(168,0.10) upper tail = 4.85e-11 (~6.47 sigma)

16 / 168 below 0.05:
  Poisson(lambda=8.4) upper tail = 0.0125
  exact Binomial(168,0.05) upper tail = 0.0105
```

The historical `47/168` value was labeled as a Poisson tail `<1e-10`; literal
Poisson arithmetic gives `1.16e-9`, whereas the exact binomial tail is
`4.85e-11`. The historical `16/168` Poisson-tail value is also not reproduced
by the stated Poisson model. These arithmetic-label corrections do not alter
the observed 168-vector itself.

The paper text also says "11 equal bins," but its figure and reported expected
counts correspond to edges
`[0,.05,.10,.20,.30,.40,.50,.60,.70,.80,.90,1.0]`. The audit uses those
figure-compatible edges.

### 10.2 Recompute the actual historical 168-vector

```bash
python scripts/run_gwtc_population_global_null.py \
  --observed tables/gwtc4_population_observed.csv \
  --p-column scan_null_p \
  --output tables/gwtc4_population_global_result.json
```

The recovered vector gives:

```text
Observed scans        = 168
chi2                   = 95.630952
analytic uniformity p = 4.06578e-16
nominal one-sided Z   = 8.0522 sigma
count p < 0.05        = 16
count p < 0.10        = 47
```

This verifies the nominal 8.05-sigma aggregate directly from the original scan
results. It is still not a correlation-corrected global significance.

### 10.3 Verify the restored selector/scanner implementation

The outer-null generator first recomputes all 168 historical scan-max
statistics from the recovered event summary and aborts unless the saved
`k_best` and `Delta-chi2` values are reproduced to tight numerical tolerance.

```bash
python scripts/generate_gwtc4_population_null_matrix.py --verify-only
```

This requires the recovered historical event-summary file:

```text
outputs_gwtc4_wct_subsets/gwtc4_subset_event_summary.csv
```

### 10.4 Generate the dependence-preserving catalog-level null

The implemented outer null permutes finite `final_spin_median` values across
event labels while keeping each event's uncertainty weight and all non-final-
spin quantities fixed. The **same 56 subset rules and 3 final-spin pairs are
rebuilt for every null catalog**, including selectors that explicitly rank or
threshold on final spin.

This is a correlation-aware event-label permutation null. It is not a full LVK
astrophysical population model, so any surviving significance should be stated
with that null model named explicitly.

Start with a validation ensemble:

```bash
python scripts/generate_gwtc4_population_null_matrix.py \
  --outer-n 200 \
  --output-matrix tables/gwtc4_population_null_matrix.csv \
  --output-observed tables/gwtc4_population_observed_outercal.csv
```

The generator writes:

- `tables/gwtc4_population_null_matrix.csv` — one correlated 168-p-value vector
  per outer-null catalog;
- `tables/gwtc4_population_observed_outercal.csv` — the observed 168 scan-max
  statistics recalibrated against the same outer-null ensemble;
- `tables/gwtc4_population_null_delta_matrix.csv` — raw scan-max Delta-chi2
  values for every null catalog;
- `tables/gwtc4_population_null_metadata.json` — null definition, settings,
  reproduction errors and resolution floors.

The per-scan outer p-values are rank calibrated, rather than obtained by
nesting hundreds of inner null scans inside every outer catalog. This keeps the
historical selector dependence while making the global calibration tractable.

### 10.5 Compute the correlation-aware population significance

Use the outer-calibrated observed vector, not the historical independently
calibrated `scan_null_p` column:

```bash
python scripts/run_gwtc_population_global_null.py \
  --observed tables/gwtc4_population_observed_outercal.csv \
  --p-column outer_global_p \
  --null-matrix tables/gwtc4_population_null_matrix.csv \
  --output tables/gwtc4_population_global_result.json
```

That final command compares the observed aggregate statistic with the aggregate
statistics of complete null-catalog 168-vectors. The reported empirical p-value
therefore includes the overlap/correlation induced by the historical selector
architecture under the stated final-spin event-label permutation null.

**Do not construct the null matrix from 168 independent Uniform(0,1) draws.**
That would erase the selector dependence and return to the independence-based
analytic approximation.

The catalog-level empirical p-value uses the plus-one correction
`(1 + N_ge) / (1 + N_null_catalogs)`. Directly resolving a one-sided 5-sigma
p-value with zero exceedances requires roughly **3.49 million complete outer
null catalogs**. Smaller ensembles still measure whether the nominal 8.05-sigma
aggregate collapses substantially once selector dependence is included.

## Outputs

- `tables/gwtc_*.csv` — canonical current-harness results (see METHOD.md §6).
- `outputs/summary/VERDICT.txt` — current diagnostic verdict.
- `outputs/summary/scan_*.csv`, `outputs/summary/*.png` — scan curves / plots.
- `tables/gwtc4_population_observed.csv` — restored historical 168-scan results.
- `tables/gwtc4_population_global_result.json` — nominal and/or correlated
  population calibration result.
- `tables/gwtc4_population_null_matrix.csv` — correlated outer-null p vectors.
- `tables/gwtc4_population_null_delta_matrix.csv` — corresponding raw scan-max
  Delta-chi2 vectors.

## Seeds

Default seeds are explicit: primary/stress `--seed 12345`, controls
`--seed 2024`, historical outer population null `--seed 20260831`. Change them
to test stability; a strong result should not hinge on one random seed.
