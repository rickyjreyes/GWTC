# GWTC Baseline-v2 and Holdout Closure Plan

This branch upgrades the current Class-III log-domain residual diagnostic without changing the historical result. The goal is to test whether the candidate residual survives stronger astrophysical baselines and a frozen holdout protocol.

## Scientific rules

1. Preserve the current published/result tables as the v1 diagnostic record.
2. Do not promote the current candidate above PARTIAL / Class III unless all v2 gates pass.
3. Freeze the primary variable, coordinate, search range, baseline family, model-selection rule, and holdout split before inspecting holdout residuals.
4. A flexible baseline must be selected using training data only.
5. The holdout test must evaluate the frozen residual frequency and phase without rescanning the holdout for a new best mode.
6. Report negative results unchanged.

## Baseline-v2 challenge

Compare the current degree-3 Poisson baseline against stronger smooth alternatives on the training set:

- polynomial degrees 3 through 8;
- cubic spline Poisson GLM bases with declared knot counts;
- Gaussian-mixture density baselines in log mass;
- cross-validated kernel-density baselines where normalization is explicit.

Baseline selection must use predictive fit on training folds. The oscillatory residual is evaluated only after the baseline family and hyperparameters are frozen.

## Unbinned cross-check

Add an unbinned likelihood diagnostic on `ell = ln(z)` so the principal result no longer depends on a histogram partition. The binned and unbinned analyses are separate checks; agreement strengthens the candidate, disagreement is reported as fragility.

## Frozen holdout

Default chronological protocol:

- development/training: earlier GWTC releases / observing runs;
- holdout: later events not used to tune the analysis;
- learn candidate `k`, phase, and amplitude on training only;
- evaluate the frozen residual on holdout without rescanning `k`;
- record predictive log-likelihood improvement and a null-calibrated holdout p-value.

If release boundaries change, record the exact event names and catalog versions in a manifest before evaluation.

## Required outputs

- `tables/gwtc_baseline_v2_cv.csv`
- `tables/gwtc_unbinned_crosscheck.csv`
- `tables/gwtc_holdout_manifest.csv`
- `tables/gwtc_holdout_result.csv`
- `outputs/summary/BASELINE_V2_VERDICT.txt`

## Promotion gates

A result may advance beyond Class III only if all of the following are true:

- candidate remains significant after the selected flexible baseline;
- inferred winding/frequency is stable under binning and unbinned analysis;
- bounded non-scale controls are not equally significant;
- the frozen mode improves prediction on held-out events;
- selection/posterior uncertainty checks do not destroy the result;
- no holdout hyperparameter or frequency was chosen after seeing holdout outcomes.

Until then, the repository remains a diagnostic harness and the current Class-III interpretation remains authoritative.
