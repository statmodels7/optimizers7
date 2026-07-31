# Limited-Memory BFGS

BFGS without ever forming the matrix: only the last `memory` secant
pairs are kept, and the direction comes from the two-loop recursion.

## Usage

``` r
lbfgs(
  criterion = crit_any(crit_grad(1e-08), crit_rel_obj(1e-12)),
  memory = 10,
  curv_tol = 1e-10,
  step = 1,
  line_search = wolfe(),
  maxit = 500,
  max_eval = 20000,
  verbose = FALSE,
  refresh = 10,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule.

- memory:

  How many secant pairs to keep. Defaults to 10.

- curv_tol:

  A pair is discarded when \\s^\top y \le \texttt{curv\\tol}\\\lVert
  s\rVert\lVert y\rVert\\. Defaults to `1e-10`.

- step, line_search, maxit, max_eval, verbose, refresh, keep_trace:

  As in
  [`bfgs`](https://statmodels7.github.io/optimizers7/reference/Bfgs.md).

## Value

An S7 object of class `Lbfgs`.

## Details

The cost is \\O(mp)\\ in both time and memory rather than \\O(p^2)\\.
For a handful of parameters that is no gain whatever and
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/Bfgs.md)
will converge in fewer iterations, since it carries the whole
approximation; the crossover is somewhere in the hundreds of parameters,
and beyond a few thousand the full matrix is simply not storable.

The recursion is scaled at each iteration by the most recent pair's
\\s^\top y / y^\top y\\. That single number is doing the work the full
matrix would otherwise do, and it is the reason the method converges at
all without one.

A larger `memory` is not uniformly better: old pairs describe curvature
at points the iterate has left, and on a strongly nonlinear objective
they can be worse than no information. Ten is the conventional choice
for good reason.

## See also

[`bfgs`](https://statmodels7.github.io/optimizers7/reference/Bfgs.md)

## Examples

``` r
rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
minimize(lbfgs(memory = 5), rosen, c(-1.2, 1))
#> <optimizer_result> L-BFGS
#>   value      : 5.34216e-17
#>   par        : 1 1
#>   iterations : 34   evaluations: f 253, g 0
#>   converged  : yes (gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative))
#>   note       : gradient obtained by finite differences
```
