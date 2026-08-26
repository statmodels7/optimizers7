# Print Method for an Optimization Result

Shows the run in six lines at most: the method's name, the objective
value, the leading parameters, the iteration and evaluation counts, the
elapsed time and the convergence status with the rule that fired. A
non-empty `message` adds a `note` line. A failure prints `NO` in
capitals, so a run that did not converge cannot be skimmed past.

## Arguments

- x:

  An
  [`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

- digits:

  Decimal places the parameters are rounded to. A single non-negative
  whole number, default 4. Anything else raises an error. The objective
  value is not affected: it always prints to six significant figures.

- max_par:

  How many parameters to show. A single positive whole number, default
  6; the remainder is reported as `... (6 of 40 shown)`.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the output.

## Examples

``` r
res <- minimize(gd(), function(p) sum((p - 1:2)^2), c(0, 0))
res
#> <optimizer_result> gradient descent
#>   value      : 2.80957e-22
#>   par        : 1 2
#>   iterations : 1   evaluations: f 11, g 0
#>   elapsed    : 1e+03 us
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
#>   note       : gradient obtained by finite differences
print(res, digits = 2, max_par = 1)
#> <optimizer_result> gradient descent
#>   value      : 2.80957e-22
#>   par        : 1 ... (1 of 2 shown)
#>   iterations : 1   evaluations: f 11, g 0
#>   elapsed    : 1e+03 us
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
#>   note       : gradient obtained by finite differences

# A run stopped by its budget says so on the converged line.
rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
print(minimize(gd(maxit = 5), rosen, c(-1.2, 1)))
#> <optimizer_result> gradient descent
#>   value      : 4.10215
#>   par        : -1.0203  1.0553
#>   iterations : 5   evaluations: f 77, g 0
#>   elapsed    : 0 us
#>   converged  : NO (iteration budget reached)
#>   note       : gradient obtained by finite differences
```
