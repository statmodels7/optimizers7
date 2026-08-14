# Summary Method for an Optimization Result

Summary Method for an Optimization Result

## Arguments

- object:

  An
  [`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

- ...:

  Unused.

## Value

`object`, invisibly. Called for the printed summary.

## Examples

``` r
res <- minimize(gd(keep_trace = TRUE),
                function(p) sum((p - 1:2)^2), c(0, 0))
summary(res)
#> <optimizer_result> gradient descent
#>   value      : 2.80957e-22
#>   par        : 1 2
#>   iterations : 1   evaluations: f 11, g 0
#>   elapsed    : 1e+03 us
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-12 (relative))
#>   note       : gradient obtained by finite differences
#>   safeguards :
#>     step shortened: 1
```
