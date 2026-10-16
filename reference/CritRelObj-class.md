# S7 Class for the Relative Objective Criterion

This is the class of the rule built by
[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md).
The rule fires when \\\lvert f\_{new} - f\_{old}\rvert \< \texttt{tol}\\
(\lvert f\_{old}\rvert + \texttt{tol})\\, so the change is compared with
the size of the objective and one tolerance serves for any units. It
declares nothing through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md).

## Usage

``` r
CritRelObj(label = character(0), tol = integer(0))
```

## Arguments

- tol:

  The tolerance, a single positive number.

## Value

An S7 object of class `CritRelObj` inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
carrying `label` and `tol`.

## Details

[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
returns `FALSE` when `f_old` is `NULL` or not finite.

The term `+ tol` is a floor on the scale. Without it the threshold would
be zero at an objective value of zero, and since the comparison is
strict the rule could never fire at an optimum of zero.

## See also

[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
for the constructor,
[`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)
for the version in the units of the objective.

## Examples

``` r
crit_met(crit_rel_obj(1e-6), list(f_new = 1.0000001, f_old = 1.0000002))
#> [1] TRUE

# With the floor, the rule can fire at an optimum of zero.
crit_met(crit_rel_obj(), list(f_new = 1e-30, f_old = 0))
#> [1] TRUE
```
