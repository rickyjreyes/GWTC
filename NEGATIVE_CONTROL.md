# NEGATIVE CONTROL DISCIPLINE

Negative controls are first-class results here, not an afterthought. They are
**retained even when negative** and are required before any PASS-grade claim.
Run via `scripts/run_gwtc_controls.py` → `tables/gwtc_control_results.csv`.

## Why control type matters

A small `global_p` means **different things** for different controls. The
verdict logic (`scripts/make_gwtc_verdict.py`) categorizes them accordingly.

### 1. Structureless controls — `smooth_resample`, `uniform_ell`

- `smooth_resample`: replace each event value with a draw from the **smooth
  fitted baseline population** in `ell`.
- `uniform_ell`: draw events **uniformly in `ell`** over the active support.

These contain **no genuine log-domain mode**. If a structureless control reaches
a `global_p` as small as the primary, the scan is **manufacturing** significance
(e.g. via a degenerate max-over-`k` statistic), and the primary is **FAIL**.

> Observed: structureless controls return `global_p ~ 0.16-0.91` — they do
> **not** beat the primary, confirming the Poisson-bootstrap null is
> well-calibrated and the scan does not fabricate significance on noise.

### 2. Robustness control — `jitter`

Perturb each event within its **asymmetric posterior uncertainty**
(`*_lower` / `*_upper`) and re-scan. Here a **small** `global_p` is a
**positive** sign: the candidate *survives* measurement error. A jitter control
is therefore **not** counted as "beating" the primary.

> Observed: jitter `global_p ~ 0.0025-0.0075` — the candidate survives jitter.

### 3. Coordinate-stress control — `bounded_chi_eff`

Run the identical scan on `|chi_eff|`, a **bounded, non-scale** coordinate. If
this is itself significant, the detected mode is **not specific to scale-like
variables**, so the primary **cannot be Class I** and is downgraded.

> Observed: `bounded_chi_eff` `global_p ~ 0.0025` — significant. The current
> primary is therefore flagged **coordinate-non-specific** and capped at
> Class III.

## Decision summary used by the verdict

| control class      | small `global_p` means | effect on primary |
|--------------------|-------------------------|-------------------|
| structureless      | scan manufactures signal| **FAIL** if it beats primary |
| jitter             | survives measurement error | robustness **positive** (not a "beat") |
| bounded non-scale  | mode not scale-specific | **cannot be Class I** |

## Class IV

Pure-null / negative-control-grade outcomes are labelled **Class IV** and are
kept in the tables. A primary that is beaten by a structureless control, or a
non-scale coordinate that drives the apparent signal, is reported as such
rather than discarded.
