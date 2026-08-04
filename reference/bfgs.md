# BFGS

Builds an approximation to the inverse Hessian from successive
gradients, so the direction is a matrix-vector product and no second
derivatives are ever required.

## Usage

``` r
bfgs(
  criterion = crit_any(crit_grad(), crit_rel_obj()),
  curv_tol = 1e-10,
  max_skip = 5,
  step = 1,
  line_search = wolfe(),
  maxit = 500,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 10,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule.

- curv_tol:

  The update is skipped when \\s^\top y \le \texttt{curv\\tol}\\\lVert
  s\rVert\lVert y\rVert\\. Defaults to `1e-10`.

- max_skip:

  Consecutive skipped updates before the approximation is reset to the
  identity. Defaults to 5.

- step, line_search, maxit, max_eval, verbose, refresh, keep_trace:

  As in
  [`newton`](https://statmodels7.github.io/optimizers7/reference/newton.md).
  The default line search is
  [`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md);
  see Details.

## Value

An S7 object of class `Bfgs`.

## Details

The default line search is the strong Wolfe one, and that is a
requirement rather than a preference. BFGS updates its approximation
from the secant pair \\(s, y)\\ with \\s = x\_{new} - x\_{old}\\ and \\y
= g\_{new} - g\_{old}\\, and the update is only meaningful when \\s^\top
y \> 0\\. Armijo backtracking can accept a step so short that the
gradient has barely moved, giving a pair with no curvature information
in it; the Wolfe curvature condition is exactly the guarantee that this
does not happen. Run on Armijo alone, BFGS loses the property that makes
it BFGS.

When the curvature condition fails anyway the update is **skipped**
rather than applied. A small \\s^\top y\\ makes \\\rho = 1/s^\top y\\
enormous and one bad step destroys the accumulated approximation; a
stale but sound matrix is better than a fresh but corrupted one. After
`max_skip` consecutive skips there is no curvature information left
worth keeping, and the matrix is reset to the identity — the method
restarts as steepest descent and rebuilds.

The first accepted pair also rescales the identity by \\s^\top y /
y^\top y\\. Without it the first quasi-Newton step is taken with a unit
Hessian, which on a badly scaled problem has entirely the wrong
magnitude and wastes a line search discovering so.

## See also

[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
for many parameters,
[`newton`](https://statmodels7.github.io/optimizers7/reference/newton.md)
when a Hessian is available.

## Examples

``` r
rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
minimize(bfgs(), rosen, c(-1.2, 1))
#> <optimizer_result> BFGS
#>   value      : 4.25987e-18
#>   par        : 1 1
#>   iterations : 33   evaluations: f 255, g 0
#>   elapsed    : 2 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-12 (relative))
#>   note       : gradient obtained by finite differences
```
