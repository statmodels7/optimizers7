# Stop When the Objective Stops Moving (Absolute)

Builds the rule \\\lvert f\_{new} - f\_{old} \rvert \< \texttt{tol}\\.
It asks nothing of the method beyond an objective, so every optimizer
accepts it, and it is one of the three terms in the gradient methods'
default rule.

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

The tolerance is in the objective's own units, so it carries the scale
of the problem: `1e-10` is a strict demand on a log-likelihood per
observation and a loose one on a summed log-likelihood of order
\\10^{5}\\.
[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
is the version that does not.

A rule reading the objective cannot tell a stalled run from a converged
one. Where that distinction matters, ask for
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md).

## See also

[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
for the scale-free version,
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
for the test of stationarity,
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
for the default rule that holds this one.

## Examples

``` r
crit_abs_obj()
#> <criterion> |df| < 1e-10
crit_abs_obj(1e-6)
#> <criterion> |df| < 1e-06

crit_met(crit_abs_obj(1e-6), list(f_new = 1.0000001, f_old = 1.0000002))
#> [1] TRUE
```
