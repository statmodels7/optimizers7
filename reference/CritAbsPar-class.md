# S7 Class for the Absolute Parameter Criterion

This is the class of the rule built by
[`crit_abs_par()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md).
The rule fires when the largest change of a coordinate, \\\max_j \lvert
x_j^{new} - x_j^{old}\rvert\\, falls below `tol`, so the tolerance is in
the units of the parameters. It declares nothing through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md).

## Usage

``` r
CritAbsPar(label = character(0), tol = integer(0))
```

## Arguments

- tol:

  The tolerance, a single positive number.

## Value

An S7 object of class `CritAbsPar` inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
carrying `label` and `tol`.

## Details

[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
returns `FALSE` when `x_old` is `NULL`. With box constraints the
comparison is on the unconstrained scale, on which the optimizer moves,
so the same tolerance has different meanings for a variance near zero
and for one near a thousand.

## See also

[`crit_abs_par()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md)
for the constructor,
[`crit_rel_par()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_par.md)
for the version relative to each coordinate.

## Examples

``` r
crit_met(crit_abs_par(1e-6),
         list(x_new = c(1, 2), x_old = c(1, 2 + 1e-9)))
#> [1] TRUE
crit_met(crit_abs_par(), list(x_new = c(1, 2), x_old = NULL))
#> [1] FALSE
```
