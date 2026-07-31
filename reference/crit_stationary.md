# Stop When the Method's Own Measure of Progress Is Small

The stopping rule for a method that has no gradient to test.

## Usage

``` r
crit_stationary(tol = 1e-08)
```

## Arguments

- tol:

  Numeric tolerance. Defaults to `1e-8`.

## Value

A
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
object.

## Details

A gradient-based method knows it has arrived because \\\nabla f\\
vanishes. None of the derivative-free methods can use that test, and for
the non-smooth problems they exist to solve it would not be the right
test even if they could: at the minimum of \\\lvert x \rvert\\ the
subgradient you happen to evaluate is \\\pm 1\\, so
[`crit_grad`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
would sit there forever while sitting exactly on the answer.

Each such method therefore reports a non-negative scalar of its own that
goes to zero as it converges, and this rule tests that. What the scalar
*is* differs, deliberately, because the natural measure differs:

- [`nelder_mead`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md):

  the diameter of the simplex, so the tolerance is on the parameter
  scale.

- [`compass`](https://statmodels7.github.io/optimizers7/reference/compass.md):

  the poll size \\\Delta\\. This is the rule with a theorem behind it:
  the limit points of a pattern search with \\\Delta \to 0\\ are Clarke
  stationary.

- [`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md):

  the optimality estimate \\\lVert p \rVert^2 + \alpha\\, which vanishes
  exactly when zero lies in the convex hull of the collected
  subgradients with no linearisation error. Note that this is *not* the
  predicted decrease, which carries a factor of the trust parameter and
  can therefore be driven to zero by that parameter shrinking rather
  than by the point becoming stationary.

The measure appears in the trace as the `stationarity` column, so a run
can be read afterwards without knowing which method produced it.

## See also

[`nelder_mead`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
[`compass`](https://statmodels7.github.io/optimizers7/reference/compass.md),
[`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md),
[`crit_grad`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)

## Examples

``` r
crit_stationary()
#> <criterion> stationarity < 1e-08
crit_stationary(1e-10)
#> <criterion> stationarity < 1e-10
```
