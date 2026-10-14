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
