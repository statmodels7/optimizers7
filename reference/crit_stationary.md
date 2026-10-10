# Stop When the Method's Own Measure of Progress Is Small

Builds the stopping rule for a method that has no gradient to test. Each
such method reports a non-negative scalar of its own that tends to zero
as the method converges, and this rule tests that scalar. It is the
default of
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md),
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md),
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) and
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md).
The gradient-based methods report no such scalar, so they reject the
rule when the run starts.

## Usage

``` r
crit_stationary(tol = 1e-08)
```

## Arguments

- tol:

  Numeric tolerance, a single positive number. Defaults to `1e-8`. The
  quantity to which the tolerance applies depends on the method; see
  Details.

## Value

An S7 object of class
[CritStationary](https://statmodels7.github.io/optimizers7/reference/CritStationary-class.md),
inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

## Details

A gradient-based method detects convergence through the vanishing of
\\\nabla f\\. The derivative-free methods cannot use that test, and on a
non-smooth problem it is not a reliable test in any case: at every point
other than the minimum of \\\lvert x \rvert\\ the gradient is \\\pm 1\\,
so
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
fires only at an iterate that lands exactly on the minimum.

The scalar differs by method:

- [`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md):

  the largest distance from the best vertex to another vertex of the
  simplex, so the tolerance is on the parameter scale.

- [`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md):

  the poll size \\\Delta\\. For a continuously differentiable objective,
  the limit points of a coordinate pattern search with \\\Delta \to 0\\
  are stationary (Torczon 1997).

- [`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md):

  the half-width of the bracket.

- [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md):

  Corana's termination measure; see
  [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md).

- [`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md):

  the optimality estimate \\\lVert p \rVert^2 + \alpha\\, which vanishes
  exactly when zero lies in the convex hull of the collected
  subgradients with no linearization error. It is not the predicted
  decrease, which contains a factor of the trust parameter and therefore
  tends to zero when that parameter shrinks, whether or not the point is
  stationary.

The measure appears in the trace as the `stationarity` column.

## See also

[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md),
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md),
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) and
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
for the five measures,
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
for the gradient-based rule.

## Examples

``` r
crit_stationary()
#> <criterion> stationarity < 1e-08
crit_stationary(1e-10)
#> <criterion> stationarity < 1e-10

# It reads state$stationarity, and an absent measure is not a small one.
c(present = crit_met(crit_stationary(1e-6), list(stationarity = 1e-9)),
  absent  = crit_met(crit_stationary(), list(stationarity = NULL)))
#> present  absent 
#>    TRUE   FALSE 

# The measure differs by method (a size of the simplex for nelder_mead(),
# the poll size for compass()), and the trace records it in the same
# column.
q <- function(p) sum((p - c(1, 2))^2)
names(minimize(nelder_mead(keep_trace = TRUE), q, c(0, 0))@trace)
#> [1] "iteration"    "value"        "stationarity" "step"         "safeguard"   
names(minimize(compass(keep_trace = TRUE), q, c(0, 0))@trace)
#> [1] "iteration"    "value"        "stationarity" "step"         "safeguard"   
```
