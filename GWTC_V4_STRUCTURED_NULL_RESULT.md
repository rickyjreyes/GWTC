# GWTC V4 structured non-periodic population-null result

## Status

**Verdict: `PASS_ALL_STRUCTURED_NULLS`**

This file records the one-shot evaluation of the previously frozen GWTC V4
structured-population challenge. V4 is a robustness/source-discrimination test,
not an independent future-catalog replication, because GWTC-5 had already been
inspected in earlier V2/V3 work.

## Frozen signal definition

The V4 evaluator used the committed frozen artifact
`tables/gwtc_v4_frozen_structured_null.json` and did not scan or retune the
signal model on the holdout.

- variable: `M_chirp`
- coordinate: `ell = log(M_chirp)`
- fixed frequency: `k = 9.7`
- fixed cosine coefficient: `a = 0.17539820818119706`
- fixed sine coefficient: `b = 0.34311596388867205`
- fixed amplitude: `0.38534801946867075`
- fixed phase: `1.0982350914698238 rad`
- fixed normalizer: `1.051624662201406`
- holdout positive finite events: `104`
- observed holdout `Delta 2logL = 9.950185579607101`

The signal-side frequency, amplitude and phase were therefore fixed before the
V4 evaluation.

## Declared non-periodic null family

Each null scenario used a training-fitted continuous broken power law plus a
truncated Gaussian peak, multiplied by a predeclared monotone selection weight

```text
(Mchirp / Mchirp_max)^gamma
```

with `gamma = 0.0, 1.5, 2.5`.

The structured null contains no periodic/log-periodic term. Each scenario was
frozen before the holdout evaluation.

## Result

| selection gamma | null >= observed | null N | empirical p | verdict |
|---:|---:|---:|---:|---|
| `0.0` | `55` | `10000` | `0.005599440055994401` | `PASS_STRUCTURED_NULL_SCENARIO` |
| `1.5` | `1` | `10000` | `0.00019998000199980003` | `PASS_STRUCTURED_NULL_SCENARIO` |
| `2.5` | `2` | `10000` | `0.00029997000299970003` | `PASS_STRUCTURED_NULL_SCENARIO` |

Worst-case p across the three predeclared scenarios:

```text
0.005599440055994401
```

Overall verdict:

```text
PASS_ALL_STRUCTURED_NULLS
```

The machine-readable result is committed as
`tables/gwtc_v4_structured_null_result.csv`.

## What this establishes

Under the declared V4 test, the exact frozen `k = 9.7`, amplitude and phase
statistic is uncommon under all three structured **non-periodic** population
scenarios. In particular, allowing broad broken-power-law population structure,
a Gaussian population peak and the three declared monotone selection-weight
scenarios did not reproduce the observed fixed-mode statistic often enough to
cross the predeclared `p < 0.05` decision threshold.

This materially strengthens the GWTC fixed-mode evidence relative to the older
binned exploratory scan. It also extends the V3 result, which had already shown
that the frozen mode survives removal of histogram binning and replacement of
the polynomial baseline by a training-only Gaussian KDE.

## What this does not establish

V4 does not replace a full LVK hierarchical population analysis. The null is a
phenomenological structured population family and uses approximate selection
weighting rather than event-posterior propagation plus injection-calibrated
selection effects.

V4 is also not an independent future-catalog replication because GWTC-5 had
already been inspected in earlier stages of the program.

The strongest remaining conventional challenge is therefore a full
posterior-aware, selection-calibrated population model with the residual mode
kept externally frozen. The strongest future evidentiary test is a genuinely
unseen catalog evaluated once against a prediction frozen before release.
