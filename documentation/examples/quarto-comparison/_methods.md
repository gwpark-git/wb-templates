## Comparison method {#comparison-method}

The worked example uses six prescribed volume fractions $\phi$, from 0 to
0.40. All data are synthetic. Define three representative reduced columns:

- $D/D_0 = 1 + 1.2\phi$, reduced diffusivity.
- $\eta/\eta_s = 1 + 2.5\phi + 2\phi^2$, relative viscosity.
- $Z = 1 + 4\phi + 10\phi^2$, illustrative osmotic compressibility factor.

These polynomial definitions establish a reproducible reference for the
documentation example. The candidate evaluates the same polynomials by
factoring $\phi$. For each quantity and concentration the signed difference
is the candidate minus the reference. The selected absolute tolerance is
$10^{-14}$ because all outputs are dimensionless and of order one. Evaluate
whether a different scale calls for absolute and relative tolerances when
you adapt the template.
