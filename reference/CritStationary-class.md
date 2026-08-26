# S7 Class for the Stationarity Criterion

The rule
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
builds, and the two methods it implements. It reads
`state$stationarity`, the non-negative measure a derivative-free method
reports in place of a gradient, and fires when it falls below `tol`. It
declares `"stationarity"` through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
so a method that reports none refuses it when the run starts.

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
returns `FALSE` when the measure is `NULL`, empty or not finite. What
the measure *is* differs by method, so the tolerance means something
different for each: the simplex diameter for
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
the poll size for
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md),
Corana's termination measure for
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md), and
the optimality estimate \\\lVert p\rVert^{2} + \alpha\\ for
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md).
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
carries the comparison across all four.

## See also

[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for the constructor and what each method reports,
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
