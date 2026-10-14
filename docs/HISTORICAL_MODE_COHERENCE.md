# Historical mode recurrence, complete-catalog outer nulls (#10)

The prior N0–N3 outer-null files store full vectors of *maximum statistic*
or calibrated p, but **not** the corresponding best k, phase, n_obs and
per-null dynamic event membership. Those matrices cannot establish mode
coherence. The recovered historical original observed CSV holds 168 scans
with k, phase, active width and event list. Freeze the new metric in
`protocols/historical_coherence_v1.json` BEFORE new outer-null scoring.

The metric deduplicates labels with identical physical-event membership,
restricts comparisons to the same (x_var,y_var), requires overlap
Jaccard>=0.5, and uses a Jaccard-weighted fraction of pairs with
|k_i-k_j|<=0.5 and circular phase distance<=pi/4. n_obs is checked against
k*active_width/(2*pi), not counted as independent evidence.

After restoring the original historical event summary into the locally
expected path, preserve the original N0-null design, selections, grid and
weights. Generate a NEW complete outer-null details artifact (do not rewrite
prior frozen files):

```bash
python scripts/generate_gwtc4_population_null_matrix.py \
 --event-summary LOCAL_RECOVERED_EVENT_SUMMARY.csv \
 --outer-n 1000 --seed 20260831 \
 --output-matrix results/coh_p.csv --output-observed results/coh_obs.csv \
 --output-delta results/coh_delta.csv --metadata results/coh_metadata.json \
 --output-coherence results/coh_full_peak_details.csv \
 --output-observed-coherence results/coh_observed_details.csv
python scripts/historical_mode_coherence.py \
 --observed results/coh_observed_details.csv \
 --null-details results/coh_full_peak_details.csv \
 --output results/coherence_result.json
```

Every outer catalog reruns the same 168 complete scan definitions; dynamic
spin selectors update their actual event membership. New detail output records
the observed-comparable peak k, phase convention atan2(-b,a), n_obs and
membership. No null realization or negative result can be dropped. This is a
historical exploratory calibration; it does not modify V1–V4 or turn the
historical nominal 8.0522 sigma into an independent discovery significance.
The original recovered event-summary input is NOT embedded in the public
repo and the 1,000-catalog run is not performed as part of this PR.

Both observed and outer-null peak exports must be requested together. Recomputed observed details are verified against the original saved k, Delta, event membership and phase. The frozen outer-null minimum depth and complete contiguous replicate IDs are enforced rather than producing an apparent definitive verdict from a short diagnostic.
