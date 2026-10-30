# Summary Method for an Optimization Result

Prints what
[`print.optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/print.optimizer_result.md)
shows, followed by a count of each safeguard that the run applied to its
own steps. The counts show how the run proceeded: a Newton run that
reports `step shortened` had steps shortened by the line search, one
that reports `hessian modified` was repairing an indefinite Hessian, and
a Barzilai-Borwein run that reports `bb curvature reset` had secant
pairs without usable curvature.

## Arguments

- object:

  An
  [`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

- ...:

  Unused.

## Value

`object`, invisibly. Called for the printed summary.

## Details

The safeguard table needs a trace, so the optimizer must have been built
with `keep_trace = TRUE`. Without a trace the method prints what
[`print.optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/print.optimizer_result.md)
prints. When every row of the trace has `safeguard` equal to `"none"`,
the output reads `safeguards : none fired`.

## Examples

``` r
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))

# On the curved valley the line search of newton() shortens a few steps.
summary(minimize(newton(keep_trace = TRUE), f, c(-1.2, 1), gr = g))
#> <optimizer_result> Newton
#>   value      : 3.74527e-21
#>   par        : 1 1
#>   iterations : 21   evaluations: f 29, g 106
#>   elapsed    : 2 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
#>   safeguards :
#>     step shortened: 4

# On a quadratic the safeguards do not fire.
summary(minimize(newton(keep_trace = TRUE),
                 function(p) sum((p - 1:2)^2), c(0, 0)))
#> <optimizer_result> Newton
#>   value      : 0
#>   par        : 1 2
#>   iterations : 2   evaluations: f 47, g 0
#>   elapsed    : 1 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
#>   note       : gradient obtained by finite differences
#>   safeguards : none fired
```
