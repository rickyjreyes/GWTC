# Reproduce the GWTC Baseline-v2 Challenge

This sequence is intentionally ordered so the holdout assignment is frozen before any holdout residual is inspected.

```bash
pip install -r requirements.txt
pytest -q

# Rebuild/fetch the canonical catalog using the existing repository pipeline.
python scripts/fetch_gwtc_catalog.py
python scripts/build_gwtc_table.py

# 1. Select a stronger smooth baseline using training-data predictive deviance only.
python scripts/run_baseline_v2_cv.py \
  --variable M_chirp --bins 40 --degrees 3,4,5,6,7,8 --folds 5 \
  --p-astro-min 0.5

# 2. Freeze the chronological train/holdout membership before residual evaluation.
python scripts/make_holdout_manifest.py \
  --holdout-prefix GWTC-5 --p-astro-min 0.5
```

The files `tables/gwtc_baseline_v2_cv.csv` and `tables/gwtc_holdout_manifest.csv` are evidence artifacts. Commit or archive their hashes before adding/running any holdout residual evaluator.

## Not yet implemented on this branch

The following remain required before the v2 scientific question is answered:

- non-polynomial flexible background families (spline/mixture/KDE);
- unbinned likelihood cross-check;
- selection-function correction;
- full posterior-sample propagation;
- frozen-frequency/phase holdout evaluator;
- null-calibrated holdout significance;
- updated `BASELINE_V2_VERDICT.txt`.

Until those stages are complete, the historical Class-III result remains the repository's scientific verdict.
