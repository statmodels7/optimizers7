# Newton's Method with a Modified Hessian

Solves \\H d = -g\\ for the direction, repairing \\H\\ when it is not
positive definite, and then line searches along it. Converges
quadratically near a minimum: on Rosenbrock from the customary start it
reaches a value of `3.7e-21` in 21 iterations, 29 objective evaluations
and 21 Hessians.

## Usage

``` r
newton(
  criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
  hessian_mod = c("eigen", "ridge"),
  floor = 1e-08,
  step = 1,
  line_search = armijo(),
  maxit = 200,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 10,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object. Defaults to
  `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`; see
  [`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
  for what that disjunction buys and costs.

- hessian_mod:

  How to repair an indefinite Hessian: `"eigen"` (default) or `"ridge"`.
  Partial matching applies and any other string is refused, naming both.

- floor:

  The smallest eigenvalue the repaired Hessian is allowed, a single
  positive number. Defaults to `1e-8`, relative to the largest
  eigenvalue.

- step:

  Initial step length offered to the line search, a single positive
  number. Defaults to `1`, which is the natural Newton step and is
  accepted unchanged near the solution.

- line_search:

  A
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object. Defaults to
  [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md),
  the curvature condition
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  needs being unnecessary here because the curvature is read from the
  Hessian.

- maxit, max_eval, verbose, refresh, keep_trace:

  As in
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).
  `maxit` defaults to 200 here where
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  uses 500, a Newton run that has not arrived in 200 iterations being in
  trouble of another kind.

## Value

An S7 object of class
[Newton](https://statmodels7.github.io/optimizers7/reference/Newton-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Why the Hessian has to be repaired

\\H^{-1}g\\ is a descent direction only when \\H\\ is positive definite.
Where the objective curves downwards the unmodified step points towards
a saddle or a maximum, and no line search can rescue it: every step
along an ascent direction increases the objective. That is the ordinary
situation far from the solution, not an edge case.

Both repairs begin with an attempted Cholesky factorization, which when
it succeeds is at once the test for positive definiteness and the solve.
When it fails:

- `"eigen"`:

  decompose \\H\\ and raise every eigenvalue below `floor` to it. The
  direction is then the Newton one in the subspace where the curvature
  is trustworthy, and gradient-like in the rest. Costs a symmetric
  eigendecomposition and gives the best-conditioned repair.

- `"ridge"`:

  add \\\tau I\\ with \\\tau\\ doubling until the factorization
  succeeds. This is Levenberg's idea: it interpolates between the Newton
  step at \\\tau = 0\\ and a scaled steepest-descent step for large
  \\\tau\\. Cheaper, and blunter.

Which repair fired is recorded in the trace, under the names
`hessian modified` and `hessian modified (capped)`, so a run that spent
its time repairing can be told from one that spent it converging.

## A repaired step is capped and an unrepaired one is not

Where the factorization succeeds, \\H^{-1}g\\ is the Newton step and
`step = 1` is its own unit, so it is taken as it stands. Where it fails,
the component of the direction in the floored subspace is
\\g_i/\lambda\_{\mathrm{floor}}\\, whose size is set by `floor` and by
no curvature of the objective. The direction is therefore scaled to
\\\min(1, 1/\lVert d\rVert\_\infty)\\, a displacement of order one in
the parameters. That is the rule
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md)
applies when its secant pair reports no curvature and the one
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
apply to their first direction, and it can only ever shorten a step, so
a run whose repaired steps were already of order one is untouched. The
same scaling reaches the gradient the method falls back on when a solve
fails.

Without the cap the length of a repaired step is unbounded and the line
search is the only thing bounding it, at one objective evaluation per
backtrack. Measured on a marginal criterion whose outer Hessian is
indefinite at ordinary points, a gradient of 29 along a direction
floored at \\10^{-8}\lambda\_{\max}\\ gave a step of \\4\times 10^{4}\\
on a log scale, and each of the 23 backtracks that followed was a whole
penalized refit that could not be evaluated.

## On the Hessian itself

When `he` is not supplied to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
the Hessian is obtained by differencing the gradient. With an analytic
gradient that is one numerical differentiation and is acceptable. When
the gradient is *also* differenced it is two composed, the one place in
the package where that happens, and the cost shows: one Newton iteration
on a quadratic in 20 unknowns takes 1682 objective evaluations without a
gradient and 42 gradient evaluations with one. An objective with neither
derivative is better served by
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
which needs no Hessian at all.

## References

Gill, P. E., Murray, W. and Wright, M. H. (1981). *Practical
Optimization*. Academic Press, London.

Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*, 2nd
edition. Springer, New York.

## See also

[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
when no Hessian is available,
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
for the line search,
[`summary.optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/summary.optimizer_result.md)
for the safeguard counts.

## Examples

``` r
rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
rosen_gr <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                          200 * (p[2] - p[1]^2))
rosen_he <- function(p) matrix(
  c(2 - 400 * (p[2] - 3 * p[1]^2), -400 * p[1],
    -400 * p[1], 200), 2, 2)

minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
#> <optimizer_result> Newton
#>   value      : 3.74398e-21
#>   par        : 1 1
#>   iterations : 21   evaluations: f 29, g 22
#>   elapsed    : 2 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)

# A saddle, where the unrepaired step points the wrong way. The Hessian at
# the start has eigenvalues 2 and -1.99, and both repairs reach the minimum
# in five iterations, reporting the repair in the trace.
sad <- function(p) p[1]^2 - p[2]^2 + 0.1 * p[2]^4
sg  <- function(p) c(2 * p[1], -2 * p[2] + 0.4 * p[2]^3)
sh  <- function(p) matrix(c(2, 0, 0, -2 + 1.2 * p[2]^2), 2, 2)
eigen(sh(c(0.5, 0.1)), only.values = TRUE)$values
#> [1]  2.000 -1.988

r <- minimize(newton(keep_trace = TRUE), sad, c(0.5, 0.1), gr = sg, he = sh)
c(value = r@value, iterations = r@iterations)
#>      value iterations 
#>       -2.5        5.0 
unique(r@trace$safeguard)
#> [1] "hessian modified (capped)" "none"                     

# summary() counts the repairs, which is how a struggling run is read.
summary(minimize(newton(keep_trace = TRUE), rosen, c(-1.2, 1),
                 gr = rosen_gr, he = rosen_he))
#> <optimizer_result> Newton
#>   value      : 3.74398e-21
#>   par        : 1 1
#>   iterations : 21   evaluations: f 29, g 22
#>   elapsed    : 3 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
#>   safeguards :
#>     step shortened: 4
```
