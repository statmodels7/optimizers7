# Stop When the Objective Stops Moving (Absolute)

Builds the rule \\\lvert f\_{new} - f\_{old} \rvert \< \texttt{tol}\\.
It reads only the objective, so every optimizer accepts it, and it is
one of the three terms of the default rule of the gradient methods.

## Usage

``` r
crit_abs_obj(tol = 1e-10)
```

## Arguments

- tol:

  Numeric tolerance, a single positive number. Defaults to `1e-10`.

## Value

An S7 object of class
[CritAbsObj](https://statmodels7.github.io/optimizers7/reference/CritAbsObj-class.md),
inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

## Details

The tolerance is in the units of the objective, so it depends on the
scale of the problem: `1e-10` is a strict condition on a log-likelihood
divided by the number of observations and a loose one on a summed
log-likelihood of order \\10^{5}\\.
[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
is the relative version.

A rule on the objective cannot distinguish a stalled run from a
converged one. Where that distinction matters, use
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md).

## See also

[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
for the relative version,
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
for the test of stationarity,
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
for the default rule that contains this one.

## Examples

``` r
crit_abs_obj()
#> <criterion> |df| < 1e-10
crit_abs_obj(1e-6)
#> <criterion> |df| < 1e-06

crit_met(crit_abs_obj(1e-6), list(f_new = 1.0000001, f_old = 1.0000002))
#> [1] TRUE
```
