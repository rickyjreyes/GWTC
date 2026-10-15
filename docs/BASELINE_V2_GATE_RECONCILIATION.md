# Baseline-v2 original acceptance reconciliation (#3)

This is an audit against the ORIGINAL issue acceptance criteria, not a rewrite
of thresholds after seeing later data.

| Original criterion | Evidence in main/this branch | Current disposition |
|---|---|---|
| CV on canonical M_chirp **training** | `scripts/run_baseline_v2_cv.py`, committed `tables/gwtc_baseline_v2_cv.csv`, strict training-only PR #5 | Implemented/committed |
| Non-polynomial smooth baseline chosen on training-only predictive fit | V3 KDE LOO bandwidth selection; `tables/gwtc_v3_frozen_unbinned_kde_mode.json` | Implemented/committed |
| Unbinned cross-check | V3 unbinned evaluator and original V3 result | Implemented/committed |
| Freeze holdout manifest before residual inspection | V1 manifest and frozen model paths; see V2/V3 reproduction records | Implemented for declared V1 protocol, with caveat earlier pre-strict CV saw the holdout |
| Evaluate training-selected mode on holdout WITHOUT rescanning k | V1 strict fixed model and holdout result `p=0.0012998700` | Implemented/committed |
| Propagate posterior/measurement uncertainty | #11 per-event provenance-verified posterior sensitivity tool added in this PR; original official HDF5 data not downloaded or run | **OPEN: scientific evidence missing** |
| Preserve PARTIAL/Class III unless *every* promotion gate passes | README/RESULTS/METHOD historical interpretation and tests/test_no_overclaim.py | Keep existing classification; do NOT infer Class I from V1/V3/V4 |
| Record negative outcomes unchanged | historical N0–N3 artifacts, V1–V4 results retained | Implemented and must remain |

The V4 structured population null (worst-case p=0.00559944) adds useful
source-discrimination robustness, but still uses simplified selection weights.
It does **not** substitute for #8 full PE/injection-calibrated hierarchical
inference or #9 genuinely unseen catalog. Do not close #3 until the original
measurement uncertainty gate is executed and verified; report a negative
posterior result as negative if obtained.

Source data note: GWTC-5.0 PE stable release version 9 is available in official
Zenodo parts 1 and 2. Retrieval, checksums and exact waveform/product choice
must be documented rather than guessed.


## Machine-readable original acceptance audit

Run `python scripts/gwtc_baseline_v2_acceptance.py --output
results/baseline_v2_gate.json` to record archived gate SHA values and a
NOT_RUN measurement-uncertainty gate. After #11 official data execution, pass
`--posterior-result PATH --posterior-sha256 EXACT_HASH` to check the complete
coverage, exact V4 freeze and original posterior-diagnostic schema. The output
always preserves `PARTIAL_CLASS_III_NOT_PROMOTED`; full scientific promotion
still requires independent provenance/interpretation rather than code status.
Never confuse a posterior-draw sensitivity diagnostic with full hierarchical
PE and injection analysis or a genuinely unseen future replication.


## Stronger archived-input consistency checks (without gate promotion)

`scripts/gwtc_baseline_v2_acceptance.py` now verifies that the selected CV mean row is degree 7 on training data, the V3 non-polynomial KDE was chosen by training-only LOO, V3 original manifest/table SHA values match, V2/V3 evaluated frequencies equal their archived frozen frequencies, the V2 B40 frequency is below its bin-edge Nyquist limit, the training/holdout counts are 173/104, the canonical holdout manifest has no duplicate names, and V3 explicitly identifies itself as robustness rather than untouched independent replication. Contradictions fail closed; missing archived sources are reported without fabricating a PASS. The scientific status remains `PARTIAL_CLASS_III_NOT_PROMOTED` and #11 remains mandatory.


The real-archive consistency gate ran in CI after explicit frozen Windows newline verification. Both V3 frozen input CSVs can be verified against original expected SHA through a **byte-only LF/CRLF transformation**, with raw checkout SHA and matching original SHA recorded separately. This does not change the event table, holdout identities or frozen source hashes. The archived V2 selected 40-bin mode at k=9.602325620315224 lies below k_N=28.273762790192603. The actual gate remains `ORIGINAL_ACCEPTANCE_NOT_CLOSED` because no official #11 posterior diagnostic was supplied.
