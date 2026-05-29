# METHOD

This document defines the WCT-motivated log-domain residual diagnostic applied
to GWTC event-catalog variables. It is a **diagnostic harness**; it does not
prove WCT and does not replace LVK population inference.

## 1. Admissible scale variables

Only **physically scale-like** variables (positive, with meaningful ratios) may
serve as **primary** evidence:

| variable   | meaning                          | source-frame |
|------------|----------------------------------|--------------|
| `M_chirp`  | chirp mass (primary recommended) | yes          |
| `M_total`  | total mass                       | yes          |
| `D_L`      | luminosity distance              | n/a          |
| `M_final`  | final mass (if available)        | yes          |
| `redshift` | redshift (if available)          | n/a          |

`E_rad` and peak luminosity are not present in the GWOSC `allevents` summary
table and are therefore reported as unavailable (INCOMPLETE for those).

Bounded ratio/spin variables (mass ratio `q`/`eta`, `chi_eff`) are **not true
scale variables**. `|chi_eff|` is included only as a **bounded-coordinate
negative control / coordinate-stress** test and can never produce a primary
PASS.

## 2. Log coordinate

For a selected variable `z > 0`:

```
ell_i = ln(z_i)
```

## 3. Event-density (histogram) construction

The retained values define the **active log-domain support**
`[ell_min, ell_max]`, with width:

```
Delta_ell_A = ell_max - ell_min
```

`ell` is binned into `B` equal-width bins giving counts `y_j` at bin centers
`ell_j`. `B` is varied in the bin-count stress test (20, 30, 40, 50).

## 4. Baseline model

A smooth Poisson baseline is fit by a log-rate GLM with a low-order polynomial
basis (default degree 3) in centered/scaled `ell`:

```
log mu0(ell) = sum_d beta_d * t^d ,   t = (ell - mean)/std
```

fitted by iteratively reweighted least squares (IRLS).

## 5. Log-cosine (log-periodic) residual model

A single log-periodic mode is added on top of the baseline as an offset:

```
log mu(ell; k, a, b, c) = log mu0(ell) + c + a cos(k ell) + b sin(k ell)
```

For fixed `k`, `(a, b, c)` are fit by Poisson GLM with `offset = log mu0`.

## 6. Scan statistic

Using the Poisson deviance `D[y || mu] = 2 * sum( y log(y/mu) - (y - mu) )`,
the improvement at wavenumber `k` is:

```
DeltaD(k) = D[y || mu0] - D[y || mu(k)]
```

Scanning `k` over a declared grid (default 0.5 .. 40.0):

```
k_star      = argmax_k DeltaD(k)
T_obs       = max_k DeltaD(k) = DeltaD(k_star)
n_star      = k_star * Delta_ell_A / (2*pi)
```

`n_star` is the **active-domain winding number** — the number of residual
oscillations across the retained log-support. A real scale-like structure
should give a **stable `n_star`** across bin counts, not merely a stable raw
`k_star`.

A **local** chi^2 p-value (2 dof) is recorded for diagnostics only. PASS is
**never** declared from the local p-value alone.

> Note on aliasing: log-periodic modes with `k` above the per-binning Nyquist
> frequency (`~pi*B/Delta_ell_A`) are aliasing-prone and not physically
> resolvable. The bin-count stress test and the negative controls are designed
> to expose such aliasing-driven peaks as bin-fragile / coordinate-non-specific.

## 7. Null model and global p-value

The primary null is a **parametric Poisson bootstrap from the fitted baseline
`mu0`**. For each of `N_null` replicates:

1. draw `y_null ~ Poisson(mu0)`
2. refit the same baseline on the replicate
3. scan the same `k`-grid
4. record `T_null = max_k DeltaD_null(k)`

The global p-value is:

```
p_global = (1 + #{ T_null >= T_obs }) / (1 + N_null)
```

A PASS-grade claim requires `N_null >= 1000`.

## 8. Negative-control discipline

See [NEGATIVE_CONTROL.md](NEGATIVE_CONTROL.md). Briefly:

- **Structureless controls** (`smooth_resample`, `uniform_ell`): if these reach
  a `p_global` as small as the primary, the scan is *manufacturing*
  significance → FAIL.
- **jitter** (perturb within posterior uncertainty): a small `p_global` here
  means the candidate *survives* measurement error (a robustness positive); it
  is not counted as a control "beating" the primary.
- **bounded `|chi_eff|`**: if this non-scale coordinate is itself significant,
  the detected mode is **not specific to scale-like variables** and the primary
  cannot be Class I.

## 9. Look-elsewhere discipline

The full search family (variables x thresholds x bin counts x observing runs)
is recorded in `tables/gwtc_master_results.csv`. Report the per-declared-primary
`p_global`; if a claim is selected after scanning multiple branches, a
family-corrected interpretation must be applied.

## 10. Verdict rules (declared a priori)

**PASS (Class I)** requires ALL of:
- declared primary variable is physically scale-like
- primary `p_global < 0.05`
- primary `N_null >= 1000`
- bin-count stress shows **stable `n_star`** (coefficient of variation < 0.15)
- threshold stress does not destroy the signal
- structureless controls do **not** beat the primary
- the bounded non-scale coordinate is **not** equally significant
- no arbitrary coordinate transform required (primary uses `ln` of a scale var)

**FAIL** if:
- primary `p_global >= 0.05`, or
- a structureless control beats the primary, or
- the signal appears only under a single arbitrary branch.

**PARTIAL** if:
- the primary is significant but fails at least one robustness criterion
  (Class III if bin/threshold/coordinate fragile, otherwise Class II).

**INCOMPLETE** if:
- data, nulls, or provenance are missing, or required stages were not run.

## Reliability classes

- **Class I**: `p_global < 0.05`; survives primary null; `n_star` stable across
  bin counts; survives threshold stress; controls do not beat primary; not
  coordinate-fragile; variable physically scale-like.
- **Class II**: significant and interesting but missing one robustness layer.
- **Class III**: fit-able local structure but coordinate/bin/threshold fragile.
- **Class IV**: null / negative-control-grade.
