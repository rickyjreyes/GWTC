# Issue #13 — E_rad and peak luminosity availability gate

Official fields may appear only in particular GWTC posterior releases or
waveform groups, and source/detector frame and units cannot be guessed.
Run the provenance-verified field audit against the same *explicit* posterior
manifest used for #11:

```bash
python scripts/audit_official_gwtc_observables.py \
 --posterior-manifest local/official_pe_manifest.csv \
 --output results/gwtc_official_field_coverage.csv
```

This emits a separate row per event and variable, exact-name presence and
SHA-256/source/dataset provenance. It does NOT silently compute
`M_total-M_final` or use unofficial catalog approximations for radiated
energy, and it does NOT treat presence as proof of compatible physical units
or conventions. All outputs initially read analysis_ready=false. Only after
official field definitions, units, source-frame handling, cross-release
consistency, secondary scan family and calibrated nulls are locked may this
issue be declared scientifically complete.


## Independent units and frame-definition gate

Run field coverage first, then fill
`manifests/secondary_official_definitions.template.json` using independently
checked official release/source documents, exact fields, release/version,
document hashes, units and source/detector frame handling. Only after reviewing
definitions set `independent_definition_reviewed=true`. Preflight command:

```bash
python scripts/gwtc_secondary_observable_preflight.py \
 --coverage results/gwtc_official_field_coverage.csv \
 --definitions local/secondary_official_definitions.json \
 --out-prefix results/secondary_preflight
```

An exact field without documented units and frame is BLOCKED; even a reviewed
field is only CONDITIONALLY ELIGIBLE. This creates no physical observation or
p-value. Next independently freeze transforms/missingness, version consistency,
complete exploratory family and dependence-preserving null before #13 results.
