# S7 Class for the Relative Parameter Criterion

The rule
[`crit_rel_par()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_par.md)
builds. It fires when every coordinate's change, divided by that
coordinate's own size, falls below `tol`, so a parameter of order
\\10^{3}\\ and one of order \\10^{-3}\\ are held to the same number of
digits. It declares nothing through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md).

## Usage

``` r
CritRelPar(label = character(0), tol = integer(0))
```

## Arguments

- tol:

  The tolerance, a single positive number, used both as the floor and as
  the threshold.

## Value

An S7 object of class `CritRelPar` inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
carrying `label` and `tol`.

## Details

The test is \\\max_j \lvert x_j^{new} - x_j^{old}\rvert / (\lvert
x_j^{old}\rvert + \texttt{tol}) \< \texttt{tol}\\, with the same floor
[CritRelObj](https://statmodels7.github.io/optimizers7/reference/CritRelObj-class.md)
uses and for the same reason: a coordinate sitting at zero would
otherwise be divided by nothing.
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
returns `FALSE` when `x_old` is `NULL`.

## See also

[`crit_rel_par()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_par.md)
for the constructor,
[`crit_abs_par()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md)
for the version in the parameters' own units.

## Examples

``` r
# The same relative change at two very different scales.
crit_met(crit_rel_par(1e-6), list(x_new = 1000, x_old = 1000.0001))
#> [1] TRUE
crit_met(crit_rel_par(1e-6), list(x_new = 1e-3, x_old = 1e-3 + 1e-10))
#> [1] TRUE
```
