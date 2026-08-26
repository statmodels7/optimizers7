# S7 Class for the Absolute Objective Criterion

The rule
[`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)
builds. It reads `state$f_new` and `state$f_old` and fires when they
differ by less than `tol`. It declares nothing through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
every optimizer evaluating the objective, so no method refuses it.

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
returns `FALSE` when `f_old` is `NULL` or not finite, which is the state
at the first iteration, so the rule cannot fire before there are two
values to compare.

The tolerance is in the objective's own units, so this rule carries the
scale of the problem with it: `1e-10` means something different for a
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

# Nothing to compare against at the first iteration.
crit_met(crit_abs_obj(), list(f_new = 1, f_old = NULL))
#> [1] FALSE
```
