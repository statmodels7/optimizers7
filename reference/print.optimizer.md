# Print Method for Optimizers

Shows an optimizer in four lines: its name, the label of its stopping
rule, its two budgets, and the settings specific to the algorithm.
`verbose`, `refresh` and `keep_trace` are not printed.

## Arguments

- x:

  An
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  object.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the output.

## Details

The settings line describes each value according to its type. A number
prints as itself; an object with a `label`, such as a
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
or a
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md),
prints its label; an object with a `name`, such as a nested optimizer,
prints its name; any other value prints as its class in angle brackets.
A setting of a new kind therefore prints a readable value and does not
stop the method.

## Examples

``` r
# The line search is an object and shows its own label.
bfgs()
#> <optimizer> BFGS
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 500, evaluations Inf
#>   settings  : step = 1, line_search = strong Wolfe (c1 = 1e-04, c2 = 0.9), curvature_tol = 1e-10, max_skip = 5, stop_decrement = 0.001

# A derivative-free method has a different rule and different settings.
nelder_mead()
#> <optimizer> nelder-mead
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 2000, evaluations Inf
#>   settings  : step = 0.1, adaptive = TRUE, max_restarts = 3, degenerate_tol = 1e-06, simplex = <NULL>

# A wrapper names the optimizer it wraps.
multistart(bfgs(), n = 4)
#> <optimizer> multistart (BFGS)
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 4, evaluations Inf
#>   settings  : optimizer = BFGS, n = 4, starts = <NULL>, spread = 1, cores = <NULL>, distinct_tol = 1e-06
```
