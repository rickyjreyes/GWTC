# Historical coherence execution obstruction: issue #10

Classification: **BLOCKED_EXTERNAL**. The essential missing input is
`outputs_gwtc4_wct_subsets/gwtc4_subset_event_summary.csv` with verified historical provenance. No complete-catalog coherence null was generated and no calibrated coherence p-value exists from this execution.

Baseline: PR #19 at `09b3aabd77dc40d98df06dee529dfb420f3591fa`; merged PR #17 at `49b9f7849c6533605acf835f151cdcddd4ce52c9` supplied the reused scanner, coherence metric and outer-null export. All remote branches were fetched. `git log --all -- outputs_gwtc4_wct_subsets/gwtc4_subset_event_summary.csv` returned no commits; a broader event-summary history search and exact file-library search did not recover the table. An unrelated fuzzy file match was an image, not historical event data.

The existing 168-row observed scan manifest is retained. It cannot reconstruct the original per-event final-spin values, uncertainties, selectors, missingness and weights needed to rebuild the overlapping and dynamic subsets. The old delta-only matrices do not contain the required k, phase and physical-event membership. They remain unchanged and are not repurposed.

Executed command:

```bash
python scripts/record_historical_coherence_preflight.py --output NEW/coherence_preflight.json
```

This invokes the existing `generate_gwtc4_population_null_matrix.py --verify-only`, which exits 1 with:

```text
Historical event-summary file not found: outputs_gwtc4_wct_subsets/gwtc4_subset_event_summary.csv
```

`results/official_2026-10-15/coherence_preflight.json` records the exact command, return code, source/config/observed-table hashes, observed columns and the frozen protocol. Zero outer catalogs were executed. No scientific acceptance test was claimed to pass or fail based on this missing input.

After restoring the provenance-verified table, first run the verifier, then execute all 1,000 declared catalogs into fresh paths:

```bash
python scripts/generate_gwtc4_population_null_matrix.py --verify-only
python scripts/generate_gwtc4_population_null_matrix.py --outer-n 1000 --seed 20260831 --output-matrix NEW/null_p.csv --output-observed NEW/observed_outercal.csv --output-delta NEW/null_delta.csv --metadata NEW/null_metadata.json --output-observed-coherence NEW/observed_coherence.csv --output-coherence NEW/null_coherence.csv
python scripts/historical_mode_coherence.py --observed NEW/observed_coherence.csv --null-details NEW/null_coherence.csv --config protocols/historical_coherence_v1.json --output NEW/coherence_result.json
```

The frozen scanner uses 4,000 k points from 0.5 to 80, degree 2, minimum 12 events and the historical subset definitions. The coherence metric retains Jaccard>=0.5, |delta k|<=0.5, circular phase distance<=pi/4, exact-membership deduplication, and full-catalog dependence. The regenerated observed k/delta/phase/n/membership must pass the existing reproduction and parser gates before calibration. Keep #10 open until that actual experiment completes; implementation and synthetic tests alone do not satisfy its observed-versus-null table criterion.
