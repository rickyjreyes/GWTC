# Baseline-v2 original acceptance reconciliation (#3)

This is an audit against the ORIGINAL issue acceptance criteria, not a rewrite
of thresholds after seeing later data.

**2026-10-16 execution verdict: PASS_CLOSED for the original execution task.**
All 104 official positive-mass holdout PE products were verified, the frozen
posterior diagnostic ran, and its numerical output was replayed locally.
The criterion review below completes #3 without promoting the historical
**PARTIAL / Class III** scientific classification. Earlier holdout exposure
remains an irreversible limitation on prospective interpretation.

| Original criterion | Evidence in main/this branch | Current disposition |
|---|---|---|
| CV on canonical M_chirp **training** | `scripts/run_baseline_v2_cv.py`, committed `tables/gwtc_baseline_v2_cv.csv`, strict training-only PR #5 | Implemented/committed |
| Non-polynomial smooth baseline chosen on training-only predictive fit | V3 KDE LOO bandwidth selection; `tables/gwtc_v3_frozen_unbinned_kde_mode.json` | Implemented/committed |
| Unbinned cross-check | V3 unbinned evaluator and original V3 result | Implemented/committed |
| Freeze holdout manifest before residual inspection | V1 manifest and frozen model paths; see V2/V3 reproduction records | Implemented for declared V1 protocol, with caveat earlier pre-strict CV saw the holdout |
| Evaluate training-selected mode on holdout WITHOUT rescanning k | V1 strict fixed model and holdout result `p=0.0012998700` | Implemented/committed |
| Propagate posterior/measurement uncertainty | #11 complete 104-event official PE execution; `results/official_2026-10-15/posterior_execution/diagnostic.json`, paired CSV and local replay | **Executed and scientifically reviewed:** positive mean score attenuated, not eliminated; no calibrated p-value |
| Preserve PARTIAL/Class III unless *every* promotion gate passes | README/RESULTS/METHOD historical interpretation and tests/test_no_overclaim.py | Keep existing classification; do NOT infer Class I from V1/V3/V4 |
| Record negative outcomes unchanged | historical N0–N3 artifacts, V1–V4 results retained | Implemented and must remain |

The V4 structured population null (worst-case p=0.00559944) adds useful
source-discrimination robustness, but still uses simplified selection weights.
It does **not** substitute for #8 full PE/injection-calibrated hierarchical
inference or #18 genuinely unseen catalog. These remain separate obligations.
The original measurement-uncertainty gate is now executed and verified. The
current closure concerns execution and honest interpretation, not an assertion
that every physical promotion gate passed or the old holdout was untouched.

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


The real-archive consistency gate ran in CI after explicit frozen Windows newline verification. Both V3 frozen input CSVs can be verified against original expected SHA through a **byte-only LF/CRLF transformation**, with raw checkout SHA and matching original SHA recorded separately. This does not change the event table, holdout identities or frozen source hashes. The archived V2 selected 40-bin mode at k=9.602325620315224 lies below k_N=28.273762790192603. The earlier pre-PE execution returned `ORIGINAL_ACCEPTANCE_NOT_CLOSED`; that historical negative gate is preserved.

## Completed posterior execution and final scientific review

The official run [38938072587](https://github.com/rickyjreyes/GWTC/actions/runs/38938072587)
verified all 104 full raw source files (55,134,074,468 bytes total), official MD5,
computed SHA-256, preferred waveform, source-frame solar-mass field and
lossless all-sample chirp projections. The frozen statistic uses 10,000
one-draw-per-event realizations, seed 20261014, without changing k or phase.

| Quantity | Result |
|---|---:|
| Point-estimate frozen statistic | 9.9501855796 |
| Exact empirical posterior mean of the additive statistic | 8.0892161628 |
| Monte Carlo posterior-draw mean | 8.1088344680 |
| Posterior-draw 5%, median, 95% | 1.9336731089, 8.1279665195, 14.2449966225 |
| Posterior draws with positive statistic | 98.53% |
| Calibrated posterior p-value | Not computed; not supplied by this diagnostic |

The mean is lower by about 18.7% using the exact empirical-posterior expectation
(18.5% using the Monte Carlo mean). The positive fixed-mode score survives this
measurement-uncertainty diagnostic with attenuation. Its posterior spread is
not a frequentist significance calibration, hierarchical likelihood or Bayes
factor. Equal-weight Kish ESS is recorded separately from unavailable
autocorrelation-adjusted chain ESS.

`results/official_2026-10-15/baseline_v2_gate_complete.json` records the exact result hash and complete
coverage. The gate deliberately requires a scientific review rather than
automatically declaring closure; this section and `POSTERIOR_EXECUTION_2026-10-16.md`
provide that review. CV, training-only non-polynomial KDE, unbinned cross-check,
frozen within-protocol evaluation, measurement propagation and unchanged
negative/classification records are verified. Earlier pre-strict CV exposure
is explicitly retained, so #3 closure cannot be cited as an unseen prediction.
N0–N3 and V1–V4 frozen artifacts are unchanged. The full hierarchical challenge
(#8) and genuinely unseen future release (#18) remain open.

The structured manual criterion review is recorded separately in
`results/official_2026-10-15/baseline_v2_original_acceptance_review.json`.
It does not change the machine gate or rewrite its earlier negative output.
