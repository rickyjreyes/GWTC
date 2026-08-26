# GWTC V3: Non-polynomial Unbinned Robustness Test

V3 asks a narrower and harder question than V2:

> Does a training-selected log-periodic residual still predict the GWTC-5
> chirp-mass shape when the smooth background is replaced by a non-polynomial
> Gaussian KDE and the comparison is performed with an unbinned likelihood?

## Interpretation limit

GWTC-5 was already inspected in V2. Therefore V3 is **not** a new independent
catalog replication, even though all V3 fitting is restricted to the training
split. It is a predeclared robustness test against two major weaknesses of the
older analysis: polynomial baseline misspecification and histogram binning.

A V3 `PASS_ROBUSTNESS_FIXED_MODE` means the frozen V3 residual predicts the
already-declared holdout significantly better than the frozen KDE baseline
under the fixed-model Monte Carlo null. It does not establish WCT.

## Mandatory order

```bash
python -m pytest -q --basetemp=.pytest_tmp

python scripts/fetch_gwtc_catalog.py
python scripts/build_gwtc_table.py

# Recreate the same canonical primary-entry split used by V2.
python scripts/make_holdout_manifest.py \
  --holdout-prefix GWTC-5 \
  --p-astro-min 0.5

sha256sum tables/gwtc_holdout_manifest.csv

# TRAIN ONLY: select KDE bandwidth by leave-one-out likelihood, then scan/freeze
# k,a,b on the declared training split. This script cannot evaluate holdout.
python scripts/fit_frozen_gwtc_unbinned_kde_mode.py \
  --variable M_chirp \
  --bandwidth-multipliers 0.5,0.75,1.0,1.25,1.5,2.0,2.5 \
  --k-min 0.5 \
  --k-max 40 \
  --n-k 120 \
  --gh-n 40

# Freeze the exact model artifact BEFORE the evaluator is run.
sha256sum tables/gwtc_v3_frozen_unbinned_kde_mode.json
git add tables/gwtc_v3_frozen_unbinned_kde_mode.json
git commit -m "Freeze GWTC V3 unbinned KDE training model"
git push origin v3-unbinned-kde

# ONE-SHOT ROBUSTNESS EVALUATION: no bandwidth selection, no k scan, no a/b refit.
python scripts/evaluate_frozen_gwtc_unbinned_kde_holdout.py \
  --null-n 10000 \
  --seed 314159

sha256sum tables/gwtc_v3_unbinned_kde_holdout_result.csv
```

## Frozen model

The training artifact stores:

- the exact manifest and clean-table SHA-256 hashes;
- the training log-values defining the Gaussian KDE mixture;
- all bandwidth candidates and their training-only leave-one-out scores;
- the selected bandwidth;
- the declared k grid;
- the full training scan;
- frozen `k`, cosine coefficient `a`, sine coefficient `b`, amplitude, phase,
  and normalization constant `Z`;
- the equivalent winding count over the training log-domain span;
- an explicit `holdout_evaluated: false` field.

The evaluator refuses to run if the manifest or clean-table hashes differ from
those used by the freeze step.

## Unbinned model

For `ell = log(M_chirp)`, the smooth baseline is

```text
f0(ell) = Gaussian KDE trained on the training events only.
```

Bandwidth is selected by leave-one-out training log likelihood. The residual
alternative is the normalized exponential tilt

```text
f1(ell) = f0(ell) exp[a cos(k ell) + b sin(k ell)] / Z.
```

The holdout statistic is

```text
Delta 2logL = 2 sum_i log[f1(ell_i)/f0(ell_i)].
```

Because the model is frozen, this simplifies to the residual tilt minus its
normalizer; no holdout density fitting is required.

The null draws the same number of events from the exact frozen KDE mixture and
recomputes the same fixed statistic. No null replicate is refit or rescanned.

## Verdict

- `PASS_ROBUSTNESS_FIXED_MODE`: positive holdout `Delta 2logL`, `p < 0.05`, and at least 1000 null replicates.
- `FAIL_ROBUSTNESS_FIXED_MODE`: the fixed residual does not meet that criterion.
- `PARTIAL_NULL_TOO_SMALL`: fewer than 1000 null replicates.

Negative outcomes are preserved unchanged.
