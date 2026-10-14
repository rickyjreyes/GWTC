# Posterior uncertainty integration (#11) — diagnostic layer

Use official GWTC PE HDF5 files, not catalog error-bar Gaussian approximations.
Create a **complete** per-event CSV manifest with columns in exact order:
`event_name,path,dataset_path,sha256,source_url,weight_path`.
Choose and document one waveform/posterior dataset per event explicitly (e.g.
`C00:Mixed/posterior_samples` when officially available). If samples have a
separate weight dataset, supply its HDF5 path; blank means uniform weights.
The tool hashes the entire original HDF5 against the declared SHA-256 and
requires source-frame chirp mass or source-frame component masses. It refuses
to substitute detector-frame chirp mass.

```bash
python scripts/gwtc_posterior_uncertainty.py \
 --posterior-manifest local/official_pe_manifest.csv \
 --output results/gwtc_posterior_diagnostic.json
```

For each Monte Carlo realization it draws **one** posterior sample per physical
event and computes the frozen V4 residual score. It reports matched-event point
estimate score, posterior score distribution, per-event ESS, missing inputs and
coverage. Do not treat 10,000 posterior samples from one event as 10,000
independent detected mergers. This diagnostic is NOT a posterior-marginal
hierarchical Bayes factor and yields NO calibrated p-value. Full PE-prior
reweighting and injection-based selection belong to issue #8.

Official sources: GWTC-5.0 PE Part 1
https://zenodo.org/records/20348005 and Part 2
https://zenodo.org/records/20348006 ; source/version and checksums must be
recorded for each event. Do not commit multi-GB PE files.


## Complete-coverage and per-event diagnostic outputs

For the actual #11 acceptance run after validating each official event/waveform,
use `--require-complete --event-table-output results/posterior_event_comparison.csv`
with a NEW `--output` name. Incomplete manifests fail before HDF5 extraction.
The CSV reports matched canonical point contribution, independently normalized
posterior-mean contribution and their difference, weighted chirp-mass median,
per-event ESS, HDF5 SHA and exact dataset/weight paths. Do not substitute this
descriptive posterior distribution for event likelihood integration, PE-prior
reweighting, injection selection or a calibrated p-value. Check source-frame
field units against the official waveform release independently.
