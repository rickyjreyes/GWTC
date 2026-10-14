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
