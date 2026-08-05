# Print Method for Optimizers

Print Method for Optimizers

## Arguments

- x:

  An
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  object.

- ...:

  Unused.

## Value

`x`, invisibly.

## Examples

``` r
print(gd())
#> <optimizer> gradient descent
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-12 (relative)
#>   budgets   : maxit 500, evaluations Inf
#>   settings  : step = 1, line_search = Armijo backtracking (c1 = 1e-04)
```
