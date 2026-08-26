# Uniform Starting Values

Draws `npar` values, coordinate by coordinate, from the starter's `min`
and `max` on the unconstrained scale. A `min` or `max` of length one is
used for every coordinate; one of length `npar` gives each its own
range, and any other length raises an error naming both lengths.

## Arguments

- starter:

  A `UniformStart` object, read for `min` and `max`.

- npar:

  The number of parameters wanted, a positive whole number.

## Value

A numeric vector of length `npar`, drawn with
[`stats::runif()`](https://rdrr.io/r/stats/Uniform.html), so
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs it.

## Examples

``` r
set.seed(1)
starting_values(start_runif(-2, 2), 4)
#> [1] -0.9379653 -0.5115044  0.2914135  1.6328312

# A range per coordinate.
set.seed(1)
starting_values(start_runif(c(-1, -10), c(1, 10)), 2)
#> [1] -0.4689827 -2.5575220

# A length that is neither 1 nor npar is a mistake, not a request.
try(starting_values(start_runif(c(-1, -2, -3)), 2))
#> Error : 'min' must have length 1 or 2, one per parameter; it has length 3.
```
