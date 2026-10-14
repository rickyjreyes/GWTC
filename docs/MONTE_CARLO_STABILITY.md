# Issue #14 — frozen numerical stability

Run the checked-in `protocols/mc_stability_v1.json` predeclared matrix for V1,
V3 and V4. It retains the exact existing model files, observable, phase,
amplitude, k, smooth baseline, source table and holdout manifest, and varies
ONLY independent random seed and Monte Carlo depth (10,000 and 50,000).
Each mode has three predeclared seeds. All results are kept, including nonpass.

```bash
python scripts/run_frozen_mc_stability.py --out-dir results/mc_stability_v1
```

The runner invokes the original frozen evaluators with separate output paths,
records source/result hashes, reconstructs exceedance counts from exact add-one
p-values, and calculates 95% Clopper-Pearson binomial intervals. It compares
classification across the COMPLETE matrix and refuses output overwrite. It
does not pick the smallest p or allow a seed to revise scientific model choices.
The large matrix is not pre-executed in this PR: report #14 COMPLETE only when
all outputs and a stability verdict have been committed and checked.
