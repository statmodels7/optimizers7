# Backtracking Line Search with the Armijo Condition

Shrink the step until the objective decreases by enough: \$\$f(x + s d)
\leq f(x) + c_1 s\\ g^\top d .\$\$

## Usage

``` r
armijo(c1 = 1e-04, shrink = 0.5, max_step = 30)
```

## Arguments

- c1:

  Sufficient-decrease constant, in \\(0, 1)\\. Defaults to `1e-4`, the
  conventional value: it demands a decrease, but only a tiny fraction of
  what the linear model predicts, so it almost never rejects a sensible
  step.

- shrink:

  Factor applied on each backtrack, in \\(0, 1)\\. Defaults to `0.5`.

- max_step:

  Maximum backtracks before the search gives up. Defaults to 30.

## Value

A
[`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
object.

## Details

The \\c_1 s\\ g^\top d\\ term is the whole point, and dropping it —
testing merely \\f\_{new} \le f\\ — is a real defect rather than a
simplification. On a quadratic with unit step the gradient update
reflects the iterate through the minimum, leaving the objective
*exactly* unchanged; the weak test accepts it, the iterate oscillates
forever, and a stopping rule watching the objective sees no change and
reports convergence at a point that is not a minimum.

Cheap: it evaluates the objective at trial points and never the
gradient. That is enough for a method that only needs to make progress,
and not enough for a quasi-Newton method, which needs
[`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md).

## See also

[`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)

## Examples

``` r
armijo()
#> <line_search> Armijo backtracking (c1 = 1e-04)
minimize(gd(line_search = armijo(shrink = 0.2)),
         function(p) sum((p - 1:2)^2), c(0, 0))
#> <optimizer_result> gradient descent
#>   value      : 2.4818e-17
#>   par        : 1 2
#>   iterations : 39   evaluations: f 239, g 0
#>   converged  : yes (gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative))
#>   note       : gradient obtained by finite differences
```
