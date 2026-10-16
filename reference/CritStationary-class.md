# S7 Class for the Stationarity Criterion

This is the class of the rule built by
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md),
which has a method for
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
and one for
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md).
The rule reads `state$stationarity`, the non-negative measure that a
method without a gradient test reports instead, and fires when it falls
below `tol`. It declares `"stationarity"` through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
so a method that does not report this measure rejects it when the run
starts.

## Usage

``` r
CritStationary(label = character(0), tol = integer(0))
```

## Arguments

- tol:

  The tolerance, a single positive number.

## Value

An S7 object of class `CritStationary` inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
carrying `label` and `tol`.

## Details

[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
returns `FALSE` when the measure is `NULL`, empty or not finite. The
measure differs by method, so the tolerance has a different meaning for
each: the largest distance from the best vertex to another vertex of the
simplex for
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
the poll size for
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md),
the half-width of the bracket for
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md),
Corana's termination measure for
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md), and
the optimality estimate \\\lVert p\rVert^{2} + \alpha\\ for
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md).

## See also

[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for the constructor and the measure of each method,
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
for the gradient-based rule.

## Examples

``` r
crit_needs(crit_stationary())
#> [1] "stationarity"
crit_met(crit_stationary(1e-6), list(stationarity = 1e-9))
#> [1] TRUE
crit_met(crit_stationary(), list(stationarity = NULL))
#> [1] FALSE
```
