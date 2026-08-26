# S7 Class for the Empty Criterion

The rule
[`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md)
builds. Its
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
returns `FALSE` at every state without reading anything, so a run
carrying it ends only when a budget runs out and reports
`converged = FALSE`. It declares nothing through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
so no optimizer refuses it.

## Usage

``` r
CritNever(label = character(0))
```

## Value

An S7 object of class `CritNever` inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
carrying `label` alone.

## Details

It carries no `tol`, having nothing to compare, and is the only
criterion class with no property beyond `label`.
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
is the one shipped method that defaults to it.

## See also

[`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md)
for the constructor,
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
for the method that uses it.

## Examples

``` r
crit_needs(crit_never())
#> character(0)
crit_met(crit_never(), list(f_new = 0, f_old = 0, gradient = c(0, 0)))
#> [1] FALSE
```
