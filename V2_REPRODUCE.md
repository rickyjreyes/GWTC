# Reproduce the strict GWTC Baseline-v2 holdout test

The order below is mandatory.  The holdout manifest is created first, all
baseline/model selection is performed on the training split only, the fitted
mode is frozen to JSON, and only then is the GWTC-5 holdout scored.

```bash
python -m pip install -r requirements.txt
python -m pytest -q --basetemp=.pytest_tmp

python scripts/fetch_gwtc_catalog.py
python scripts/build_gwtc_table.py

# 1. Freeze the canonical cumulative-primary split BEFORE model selection.
python scripts/make_holdout_manifest.py \
  --holdout-prefix GWTC-5 --p-astro-min 0.5

# 2. Select the smooth polynomial baseline on TRAIN events only.
python scripts/run_baseline_v2_cv.py \
  --manifest tables/gwtc_holdout_manifest.csv --split train \
  --variable M_chirp --bins 40 --degrees 3,4,5,6,7,8 --folds 5 \
  --p-astro-min 0.5

# 3. Fit/scan on TRAIN only and freeze k, phase, amplitude, support, baseline.
python scripts/fit_frozen_gwtc_mode.py \
  --variable M_chirp --bins 40 --k-min 0.5 --k-max 40 --n-k 120

# Record the pre-holdout evidence hashes before evaluation.
sha256sum tables/gwtc_holdout_manifest.csv
sha256sum tables/gwtc_baseline_v2_cv.csv
sha256sum outputs/summary/gwtc_frozen_mode.json

# 4. ONE-SHOT holdout evaluation. No k scan or parameter refit occurs here.
python scripts/evaluate_frozen_gwtc_holdout.py \
  --null-n 10000 --seed 271828
```

Evidence artifacts:

- `tables/gwtc_holdout_manifest.csv`
- `tables/gwtc_baseline_v2_cv.csv`
- `outputs/summary/gwtc_frozen_mode.json`
- `tables/gwtc_frozen_holdout_result.csv`

## Interpretation

`PASS_FIXED_MODE` means only that the training-frozen residual shape predicts
the holdout counts better than the training-frozen smooth baseline under the
fixed-model null. It is **not** a WCT Class-I result.

`FAIL` means the frozen mode does not improve holdout prediction.
`PARTIAL` means the frozen mode improves prediction but does not reach the
predeclared holdout significance threshold. `INCOMPLETE_OUT_OF_SUPPORT` means
one or more positive holdout values fall outside the training-frozen log-domain
support, so no strong predictive verdict is allowed.

## Why the earlier degree-8 result must be rerun

The first baseline-v2 command used all `p_astro >= 0.5` rows before the GWTC-5
split was applied.  That was a useful baseline stress test, but it leaked
holdout events into baseline selection and therefore cannot support a strict
held-out prediction.  The sequence above repairs that before any frozen-mode
holdout result is inspected.

## Still required for a stronger scientific claim

Even a successful fixed-mode holdout remains below Class I until additional
robustness layers are passed, including a non-polynomial flexible background,
unbinned likelihood cross-check, selection-function correction, and fuller
posterior-sample propagation.
