# Print Method for an Optimisation Result

Prints the objective value, the leading parameters, the evaluation
counts, the elapsed time and the convergence status.

## Arguments

- x:

  An
  [`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

- digits:

  Decimal places the parameters are rounded to. Defaults to 4.

- max_par:

  How many parameters to show; any remainder is summarised as a count.
  Defaults to 6.

- ...:

  Unused.

## Value

`x`, invisibly.

## Examples

``` r
res <- minimize(gd(), function(p) sum((p - 1:2)^2), c(0, 0))
print(res)
#> <optimizer_result> gradient descent
#>   value      : 2.80957e-22
#>   par        : 1 2
#>   iterations : 1   evaluations: f 11, g 0
#>   elapsed    : 1 ms
#>   converged  : yes (gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative))
#>   note       : gradient obtained by finite differences
print(res, digits = 2, max_par = 1)
#> <optimizer_result> gradient descent
#>   value      : 2.80957e-22
#>   par        : 1 ... (1 of 2 shown)
#>   iterations : 1   evaluations: f 11, g 0
#>   elapsed    : 1 ms
#>   converged  : yes (gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative))
#>   note       : gradient obtained by finite differences
```
