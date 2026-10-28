# What Simulated Annealing Can Offer a Stopping Rule

Reports `"stationarity"` only. The method draws proposals and compares
values; it computes no derivative, so a rule that reads a gradient is
rejected when the run starts. Without this check the rule would test
`NULL` at every iteration and never fire.

## Arguments

- optimizer:

  An `Sa` object.

## Value

The character vector `"stationarity"`.

## Details

The measure is Corana's termination rule: at the end of each temperature
level, the largest absolute difference between the value at which the
walk stopped and both the best value found and the values at which it
stopped at the last `stationarity_levels` levels. It is read by
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md),
which is the default rule of this method.

## Examples

``` r
optimizer_provides(sa())
#> [1] "stationarity"

# A gradient rule is therefore rejected, and the message names both.
rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
try(minimize(sa(criterion = crit_grad()), rastrigin, c(1, 1)))
#> Error : The stopping rule needs gradient, which simulated annealing (uniform) does not provide.
#>   Choose a criterion that this optimizer can evaluate, or a method that provides it.
```
