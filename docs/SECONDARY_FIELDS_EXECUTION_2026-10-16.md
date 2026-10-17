# Official secondary-field availability: issue #13

**PARTIAL_OPEN.** The official quantities exist in every one of the 104 verified
preferred GWTC-5.0 posterior datasets: `radiated_energy` and `peak_luminosity`
are present by exact field name. The earlier absence-of-validated-data
obstruction has therefore narrowed. This is not a completed secondary analysis.

PR #23 builds on PR #19 at `09b3aabd77dc40d98df06dee529dfb420f3591fa` and the
merged PR #17 field gate. The full source files were size/MD5/SHA-256 verified
in workflow 38938072587. The field audit uses saved **original raw dataset**
inventories, bound to exact source identities and verified chirp projections;
it does not infer secondary fields from the derivative chirp-only datasets.
Source URLs, stable release 9 / Zenodo record identities, exact preferred
waveforms and full-file hashes are preserved in the posterior source receipts.

`secondary_coverage_complete.csv` reports 208 event/observable rows, all with
`PRESENT_EXACT_NAME`, and all with `analysis_ready=false`. The preliminary
59-event audit and failed transfer receipts are retained separately; missing
raw retrieval is distinguished from an absent field.

## Definitions found, validation still needed

The official GWOSC release documentation links the standard PESummary names.
The [parameter definitions](https://lscsoft.docs.ligo.org/pesummary/stable/gw/parameters.html)
describe radiated energy through source total and remnant mass and distinguish
ISCO-evolved from non-evolved remnant quantities. The
[PESummary unit labels](https://lscsoft.docs.ligo.org/pesummary/stable/core/file_formats.html)
give radiated energy in solar-mass equivalent and peak luminosity in
10^56 erg/s. These are documented candidate conventions, not units guessed
from observed numerical magnitudes. The two retrieved definition pages are
preserved in `secondary_definition_sources.tar.gz`, with hashes in
`secondary_definitions_review.json`.

The remaining gates are precise:

1. Bind the official release's actual producer version, remnant-fit settings
   and frame conventions to the definitions, and independently review that
   release-specific binding. Generic current documentation is not silently
   treated as a producer-version attestation.
2. Reingest and preserve the secondary values themselves, validate missing,
   nonfinite and nonpositive behavior, and check cross-version compatibility.
   The bounded primary pipeline retained only source-chirp samples; the
   presence audit is not a numerical secondary-value validation.
3. Freeze the secondary transformations, selections, frequency/phase or scan
   choices, null generator and complete exploratory-family multiplicity
   treatment **before** secondary scoring. Keep them secondary to chirp mass.

The existing fail-closed preflight was executed on all 208 rows. It reports
208 blocked, zero conditionally accepted, and no scans or null calibration.
This is an unmet analysis gate, not evidence that either physical observable
supports or refutes the residual. No mass-difference substitute was generated
and presented as an official observation. No secondary significance is claimed.

## Reproduction

Using the verified source receipts/projections described in
`POSTERIOR_EXECUTION_2026-10-16.md`:

```bash
python scripts/gwtc_secondary_receipt_audit.py --receipts NEW/pe_products --output NEW/secondary_coverage.csv
python scripts/gwtc_secondary_observable_preflight.py --coverage NEW/secondary_coverage.csv --definitions results/official_2026-10-15/secondary_definitions_review.json --out-prefix NEW/secondary_preflight
```

The complete coverage/preflight CSVs, JSON summary, definition review and
overall execution receipt are under `results/official_2026-10-15`. Existing
exact-field/unit/definition tests and the full 178-test suite pass; these tests
do not establish cross-version measurement consistency or a scientific result.
Keep #13 open until those remaining gates and the actual diagnostic execute.
