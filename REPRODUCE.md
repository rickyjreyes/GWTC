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

## 9. One-shot

```bash
bash run_all.sh    # convenience wrapper for steps 5-8 (see repo root)
```

## Outputs

- `tables/gwtc_*.csv` — canonical results (one schema; see METHOD.md §6).
- `outputs/summary/VERDICT.txt` — human-readable verdict.
- `outputs/summary/scan_*.csv`, `outputs/summary/*.png` — scan curves / plots.

## Seeds

Default seeds are explicit: primary/stress `--seed 12345`, controls
`--seed 2024`. Change them to test stability; the verdict should not hinge on a
single seed.
