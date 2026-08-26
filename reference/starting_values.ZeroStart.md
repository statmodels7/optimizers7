# Zero Starting Values

Returns `npar` zeros. Read on the unconstrained scale, so
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
maps them through the bounds before any method sees them: a positive
parameter starts at 1, a probability at 0.5, an unbounded one at 0.

## Arguments

- starter:

  A `ZeroStart` object. Its own `npar` is not read here;
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
  has already settled the count and passes it.

- npar:

  The number of parameters wanted, a positive whole number.

## Value

A numeric vector of `npar` zeros.

## Examples

``` r
starting_values(start_zeros(), 4)
#> [1] 0 0 0 0
```
