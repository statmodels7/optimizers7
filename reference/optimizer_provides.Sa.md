# What Simulated Annealing Can Offer a Stopping Rule

Reports `"stationarity"` and nothing else. The method draws proposals
and compares values; it computes no derivative, so a rule reading a
gradient is refused when the run starts instead of sitting there testing
`NULL` and never firing.

## Arguments

- optimizer:

  An `Sa` object.

## Value

The character vector `"stationarity"`.

## Details

The measure offered is Corana's termination rule: by how much the best
value has moved over the last `n_eps` temperature levels. It is read by
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md),
which is this method's default rule.

## Examples

``` r
optimizer_provides(sa())
#> [1] "stationarity"

# So a gradient rule is refused, with both names in the message.
rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
try(minimize(sa(criterion = crit_grad()), rastrigin, c(1, 1)))
#> Error : The stopping rule needs gradient, which simulated annealing (uniform) does not provide.
#>   Choose a criterion this optimizer can evaluate, or a method that provides it.
```
