# S7 Class for the Absolute Objective Criterion

The class of the rule built by
[`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md).
The rule reads `state$f_new` and `state$f_old` and fires when they
differ by less than `tol`. It declares nothing through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
because every optimizer evaluates the objective, so every method accepts
it.

## Usage

``` r
CritAbsObj(label = character(0), tol = integer(0))
```

## Arguments

- tol:

  The tolerance, a single positive number.

## Value

An S7 object of class `CritAbsObj` inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
carrying `label` and `tol`.

## Details

[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
returns `FALSE` when `f_old` is `NULL` or not finite, so the rule cannot
fire without two values to compare.

The tolerance is in the units of the objective, so the rule depends on
the scale of the problem: `1e-10` has a different meaning for a
log-likelihood of order one and for one of order \\10^{5}\\.

## See also

[`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)
for the constructor,
[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
for the scale-free version.

## Examples

``` r
crit_needs(crit_abs_obj())
#> character(0)
crit_met(crit_abs_obj(1e-6),
         list(f_new = 1.0000001, f_old = 1.0000002))
#> [1] TRUE

# Without a previous value there is nothing to compare.
crit_met(crit_abs_obj(), list(f_new = 1, f_old = NULL))
#> [1] FALSE
```
