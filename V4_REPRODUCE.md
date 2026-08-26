# GWTC V4: fixed-k structured astrophysical-null challenge

V4 asks a narrower and harder question than V2/V3:

> How often can a non-periodic structured chirp-mass population, combined with
> a monotone detection-selection weighting, produce a fixed `k = 9.7` residual
> statistic at least as large as the observed GWTC-5 statistic?

The prospective `k = 9.7` prediction was already externally published before
V4.  **V4 is not allowed to scan or retune k.**  GWTC-5 has also already been
inspected in V2/V3, so V4 is a robustness challenge, not an independent
future-catalog replication.

## Declared structured null

The intrinsic chirp-mass density is a continuous broken power law plus a
truncated Gaussian peak.  The observed density is weighted by

```text
selection(Mchirp) proportional to (Mchirp / Mchirp_max)^gamma
```

for the predeclared sensitivity scenarios

```text
gamma = 0.0, 1.5, 2.5
```

`gamma = 2.5` is an approximate low-redshift inspiral-volume stress test.  It
is **not** a full detector injection campaign or a substitute for hierarchical
LVK population inference.

## Mandatory order

### 1. Update the V4 branch and run all tests

```bash
git switch v4-astrophysical-null
git pull origin v4-astrophysical-null
python -m pytest -q --basetemp=.pytest_tmp
```

### 2. Freeze the V4 training-only challenge

Do not run the V4 evaluator first.

```bash
python scripts/freeze_gwtc_v4_structured_null.py
```

The freeze script:

- verifies the V3 clean-table and manifest hashes;
- uses only manifest rows marked `train`;
- fixes `k = 9.7` with no frequency scan;
- fits only `a,b` at that fixed frequency against the already frozen V3 KDE;
- fits each structured non-periodic population scenario on training masses only;
- writes `tables/gwtc_v4_frozen_structured_null.json`;
- does not evaluate the holdout.

### 3. Hash and commit the frozen V4 definition before holdout evaluation

```bash
sha256sum tables/gwtc_v4_frozen_structured_null.json

git add tables/gwtc_v4_frozen_structured_null.json
git commit -m "Freeze GWTC V4 structured astrophysical-null challenge"
git push origin v4-astrophysical-null
```

Record the SHA-256 value externally before proceeding.

### 4. Run the one-shot V4 robustness challenge

Only after the frozen JSON is committed:

```bash
python scripts/evaluate_gwtc_v4_structured_null.py \
  --null-n 10000 \
  --seed 2718281
```

The evaluator cannot scan `k`, refit `a,b`, refit a structured population, or
change a selection scenario.  It evaluates the actual GWTC-5 holdout once and
simulates fixed-size catalogs from each frozen structured-null scenario.

### 5. Hash and preserve the result

```bash
sha256sum tables/gwtc_v4_structured_null_result.csv

git add tables/gwtc_v4_structured_null_result.csv
git commit -m "Record GWTC V4 structured-null robustness result"
git push origin v4-astrophysical-null
```

## Decision rule

For each structured-null scenario,

```text
p_structured = (1 + number(null statistic >= observed statistic)) /
               (1 + number of null simulations)
```

The aggregate verdict is

```text
PASS_ALL_STRUCTURED_NULLS
```

only if the fixed observed statistic is positive and **every** predeclared
scenario has `p_structured < 0.05` with at least 1000 null simulations.

A failure of any scenario is retained as

```text
FAIL_AT_LEAST_ONE_STRUCTURED_NULL
```

and is not rescued by changing `k`, the selection exponents, or the population
family after looking at the result.

## Interpretation

A V4 pass would show that the already frozen `k = 9.7` residual is not easily
reproduced by this declared family of structured non-periodic mass populations
under three simple selection-weight scenarios.

It would **not** establish that selection effects are fully excluded.  The
next stronger challenge would require event posterior propagation and a full
selection-corrected hierarchical population model using detector injections or
an equivalent calibrated detection-efficiency calculation.
