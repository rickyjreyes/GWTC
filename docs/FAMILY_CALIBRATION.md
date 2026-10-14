# Full-family exploratory calibration (#16)

This is for **new exploratory scans**, not retroactive reclassification of V1–V4.
Before inspecting observed residuals, commit a manifest such as:

\`\`\`json
{"analysis_class":"exploratory","frozen_before_null":true,
 "cell_ids":["M_chirp_B20_pa05","M_chirp_B30_pa05","M_total_B20_pa05"],
 "source_release":"<official release/version>","cuts":"<exact definition>",
 "baseline":"<declared family>","k_grid":"<per-cell Nyquist-truncated grid>",
 "null_generator":"<complete-catalog procedure and seed schedule>"}
\`\`\`

Generate one null matrix with **one entire catalog** per row. Reapply ALL
selectors, backgrounds, cut branches and the matching resolved frequency grid
for each row. Columns are the declared \`cell_ids\`, storing the scan-max
statistic. Overlap/dependence lives in the rows; per-cell independent bootstrap
files are NOT a valid input.

\`\`\`bash
python scripts/calibrate_exploratory_family.py \
 --manifest manifests/exploratory.json \
 --observed observed_cells.csv --null-matrix complete_catalog_nulls.csv \
 --out-prefix results/exploratory
\`\`\`

The outputs preserve raw statistic, local empirical p, *joint family-global*
p, complete null min-p distribution, and source hashes. No significant local
cell alone is eligible for promotion to a frozen prediction. Do not optimize
the manifest or null family after observed evaluation. Retain negative cells,
null outcomes, and the provenance of every complete null realization.

For historical GWTC-4, N0–N3 already used complete-catalog calibrations; keep
those frozen historical artifacts and their explicit null sensitivity distinct.
