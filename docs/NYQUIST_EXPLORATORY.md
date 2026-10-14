# Exploratory binned Nyquist rule (issue #15)

For B uniformly spaced centers across log-domain width Delta ell, the bin
spacing is Delta ell / B and the angular Nyquist bound is
`k_N = pi B / Delta ell`. The exact boundary is excluded because sampled
sine/cosine basis columns can be degenerate there. The shared
`resolved_k_grid(field, declared_grid)` retains *only existing* frequencies
strictly below `k_N*(1-1e-10)`; it never constructs an unreported replacement
grid. Both observed and bootstrap scans consume that same filtered vector.
Direct calls to `scan_k` or `run_poisson_null` reject aliased grids.

The legacy summary notes record per-cell k_N and the number of excluded points.
Primary per-k curves record k_N. Historical outputs and original high-k scan
maxima are NOT rewritten or retroactively treated as resolved measurements.
The restored 168-scan GWTC-4 workflow has its own historical frequency grid;
historical results remain descriptive and should be annotated if reused.

This numerical sampling limit is not an astrophysical upper-frequency bound.
It does not apply to unbinned likelihoods.

## Archived read-only presentation overlays

Run `python scripts/annotate_legacy_nyquist.py --input tables/gwtc_scan_primary.csv --output results/legacy_nyquist/primary_new.csv` (choose a new output filename; no overwrite). The repository also includes frozen-base overlays at `results/legacy_nyquist/gwtc_scan_primary_annotated.csv` and `gwtc_bin_stress_annotated.csv`. The legacy width is inferred from saved `n_star=k_best*Delta_ell/(2*pi)`, not from newly retrieved events.

**Historical primary B=20 k=32.6975, n_star=23.1292 is ABOVE its binned Nyquist limit; it is NOT a resolved mode.** The B=20 bin-stress row is likewise aliased, while B=30/40/50 k=4.48319 rows are below their respective bounds. These overlays retain original reported p/verdict as historical artifacts, but those labels cannot justify a new physical inference from aliased peaks. The V3 unbinned result is not subject to binned Nyquist classification.

The original tables, figure files and p-values are not rewritten. Regenerate any publication use of these legacy binned peaks with the explicit annotation; do not display an aliased peak without this warning.
