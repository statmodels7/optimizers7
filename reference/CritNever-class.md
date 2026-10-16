# S7 Class for the Empty Criterion

This is the class of the rule built by
[`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md).
Its
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
returns `FALSE` at every state without reading it, so a run with this
rule ends when a budget runs out or when the method stops for a reason
of its own. It declares nothing through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
so every optimizer accepts it.

## Usage

``` r
CritNever(label = character(0))
```

## Value

An S7 object of class `CritNever` inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
carrying only `label`.

## Details

It has no `tol`, since it compares nothing, and it is the only criterion
class with no property beyond `label`.
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
is the only shipped method that uses it by default.

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
