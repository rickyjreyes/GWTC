# Issue #8 — full hierarchical challenge: numerical core and external gates

A correct population comparison needs the entire source-parameter vector
theta, not a KDE of chirp mass alone. The numerical core implements a
**fixed-population-hyperparameter, rate-marginalized likelihood ratio** at
the externally frozen V4 k=9.7, a_cos, b_sin. It uses PE-posterior reference
prior reweighting and recovered-injection proposal reweighting. The mode's
unknown intrinsic normalization cancels between the per-event integrals and
selection factor when the same p0(theta|Lambda) is used throughout.

Inputs are **prepared, validated importance ratios**, not raw h5 datasets:
per-event NPZ `ell, w0_over_pe_prior`, where weight is reference PE sample
weight times p0(theta|Lambda)/pi_PE(theta); recovered-injection NPZ
`ell_detected, w0_over_injection_proposal_detected`, where weight incorporates
p0(theta|Lambda)/q_inj(theta), correct campaign mixture and observing exposure.
The JSON manifest explicitly lists each NPZ original SHA-256, official source,
population model/hyperparameters, reference priors, proposal/exposure and
detection rule, and total number of generated injections. Do not approximate
the proposal by a chirp-mass marginal or skip uncertain redshift, mass-ratio,
spin and waveform dependence.

The core computes exact fixed-theta importance-sampling 2logLR, selection
ESS, PE ESS, a correlated delta-method estimate of numerical logLR variance,
and refuses numerical PASS if requirements fail. A result does **not** imply
physical verification of externally supplied ratios and reports NO p-value,
population-parameter marginal evidence, or Bayes factor. Remaining #8 gate:
- obtain and hash GWTC-5 official PE and cumulative injection data;
- independently validate model, PE prior and selection proposal densities;
- fit/freeze physically flexible no-periodic p0 on training only;
- vary/marginalize baseline population hyperparameters consistently, validate
  injection effective sample and variance, and run complete physical nulls;
- report likelihood/model comparison with uncertainties and preserve negative
  outcomes without rescanning the frozen mode.

Official data: https://zenodo.org/records/20348005 ;
https://zenodo.org/records/20348006 ;
https://zenodo.org/records/19500052 .
