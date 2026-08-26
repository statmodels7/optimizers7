# Summary Method for an Optimization Result

Everything
[`print.optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/print.optimizer_result.md)
shows, followed by a count of each safeguard the run applied to its own
steps. This is the view that answers *why* a run behaved as it did: a
Newton run that shortened four steps was repairing an indefinite
Hessian, and a Barzilai-Borwein run resetting its curvature had secant
pairs carrying none.

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
with `keep_trace = TRUE`. Without one the method prints what
[`print.optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/print.optimizer_result.md)
prints and stops there. With a trace in which nothing fired it says
`safeguards : none fired`, so the absence is reported explicitly.

## Examples

``` r
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))

# Newton on the curved valley shortens a few steps on the way in.
summary(minimize(newton(keep_trace = TRUE), f, c(-1.2, 1), gr = g))
#> <optimizer_result> Newton
#>   value      : 3.74523e-21
#>   par        : 1 1
#>   iterations : 21   evaluations: f 29, g 106
#>   elapsed    : 2 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
#>   safeguards :
#>     step shortened: 4

# A quadratic gives it no trouble at all.
summary(minimize(newton(keep_trace = TRUE),
                 function(p) sum((p - 1:2)^2), c(0, 0)))
#> <optimizer_result> Newton
#>   value      : 0
#>   par        : 1 2
#>   iterations : 2   evaluations: f 47, g 0
#>   elapsed    : 1e+03 us
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
#>   note       : gradient obtained by finite differences
#>   safeguards : none fired
```
