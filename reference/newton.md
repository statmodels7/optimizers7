# Newton's Method with a Modified Hessian

Solves \\H d = -g\\ for the direction, repairing \\H\\ when it is not
positive definite, and then performs a line search along it. The method
converges quadratically near a minimum.

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
  keep_trace = FALSE,
  max_length = Inf,
  typical = NULL
)
```

## Arguments

- criterion:

  The stopping rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object. Defaults to
  `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`; see
  [`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
  for the effect of that disjunction.

- hessian_mod:

  How to repair an indefinite Hessian: `"eigen"` (default) or `"ridge"`.
  Partial matching applies, and any other string is rejected by
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

- floor:

  The lower bound of the repair, a single positive number. Defaults to
  `1e-8`. With `"eigen"`, eigenvalues below \\\max(\texttt{floor},
  \texttt{floor}\\\max_i\lvert\lambda_i\rvert)\\ are raised to that
  value, so the bound is relative to the largest eigenvalue when that
  eigenvalue exceeds one in absolute value; with `"ridge"`, `floor` is
  the smallest shift \\\tau\\ tried.

- step:

  Initial step length offered to the line search, a single positive
  number. Defaults to `1`, which is the Newton step itself and is
  accepted unchanged near the solution.

- line_search:

  A
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object. Defaults to
  [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md).
  The curvature condition that
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  needs is unnecessary here, because the curvature is read from the
  Hessian.

- maxit, max_eval, verbose, refresh, keep_trace:

  As in
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).
  `maxit` defaults to 200 here and to 500 in
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  and
  [`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md):
  a Newton run that has not converged in 200 iterations usually has a
  problem that more iterations do not solve.

- max_length:

  The largest component that a step may have, a single positive number.
  A direction longer than this in the infinity norm is scaled down to
  this length before the line search. Defaults to `Inf`, which leaves
  every step as it is. See below.

- typical:

  `NULL` (the default), or the typical size of each parameter: one
  positive number, or one per parameter. Every length used by the method
  is then measured in these units, as the largest \\\lvert d_i\rvert /
  t_i\\ in place of the largest \\\lvert d_i\rvert\\: the bound
  `max_length`, and the scaling of a direction whose length is not set
  by the curvature. With `NULL` the lengths are measured in the units of
  the parameters.

## Value

An S7 object of class
[Newton](https://statmodels7.github.io/optimizers7/reference/Newton-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Why the Hessian has to be repaired

\\-H^{-1}g\\ is guaranteed to be a descent direction only when \\H\\ is
positive definite. Where the objective curves downwards the unmodified
step can point towards a saddle or a maximum, and a line search along an
ascent direction cannot decrease the objective. Indefinite Hessians are
common far from the solution.

Both repairs begin with an attempted Cholesky factorization, which, when
it succeeds, is both the test for positive definiteness and the solve.
When it fails:

- `"eigen"`:

  decompose \\H\\ and raise every eigenvalue below
  \\\max(\texttt{floor}, \texttt{floor}\\\max_i\lvert\lambda_i\rvert)\\
  to that value. The direction is then the Newton direction in the
  subspace where the curvature is positive, and close to the gradient
  direction in the rest. The repair costs a symmetric eigendecomposition
  and gives the best-conditioned matrix. If the eigendecomposition
  fails, the ridge repair is used.

- `"ridge"`:

  add \\\tau I\\, with \\\tau\\ doubled until the factorization
  succeeds. This is Levenberg's idea: it interpolates between the Newton
  step at \\\tau = 0\\ and a scaled steepest-descent step for large
  \\\tau\\. It costs less than the eigendecomposition and shifts every
  eigenvalue, not only the small ones.

A repair of the Hessian is recorded in the trace as `hessian modified`,
or `hessian modified (capped)` when the repaired step is also scaled
down (see below). The label does not record which of the two repairs was
used.

## Bounding the length of a step

The Newton step minimizes the local quadratic model of the objective,
and its length is set by the curvature at the current point. Where the
objective flattens towards an asymptote the curvature is small and the
step is long, and it can carry the iterate past the minimum into the
flat region. The line search accepts such a step whenever it lowers the
objective, and in the flat region the gradient is close to zero, so the
run can stop there. Such steps occur, for example, when a smoothing
parameter is optimized on the logarithmic scale.

`max_length` bounds the step; it is the `maxstep` of the line-search
Newton method of Dennis and Schnabel (1983) and the `maxNstep` of the
smoothing-parameter iteration of mgcv. The bound only shortens a step,
and near a minimum the Newton steps are short, so the quadratic
convergence of the method is kept. A step shortened by the bound is
recorded in the trace as `newton step capped`, except at an iteration at
which the Hessian was repaired, which keeps the label of the repair.

## Parameters in different units

A length in the units of the parameters has different meanings for
parameters on different scales: a displacement of one is a long step for
a parameter on the log scale and a short one for a variance of several
hundred. `typical` gives the typical size \\t_i\\ of each parameter, as
in the scaled norm of Dennis and Schnabel (1983): every length is then
the largest \\\lvert d_i\rvert / t_i\\, for the bound `max_length` and
for the scaling of a repaired direction or of a gradient fallback.

## The length of a repaired step

Where the factorization succeeds, \\-H^{-1}g\\ is the Newton step and
`step = 1` is its natural unit, so the step is taken as it stands,
subject only to `max_length`. Where it fails, the component of the
direction in the floored subspace is \\g_i/\lambda\_{\mathrm{floor}}\\,
whose size is set by `floor` and not by the curvature of the objective.
The direction is therefore multiplied by \\\min(1, 1/\lVert
d\rVert\_\infty)\\, which gives a displacement of order one in the
parameters.
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
apply the same rule to their first direction. The factor never exceeds
one, so a run whose repaired steps were already of order one is
unchanged. The same scaling applies to the gradient direction that the
method uses when a solve fails.

Without this factor the length of a repaired step is bounded only by the
line search, at one objective evaluation per backtrack. Each backtrack
can be expensive, for example a complete penalized refit when the
objective is a marginal criterion.

## The Hessian by differences

When `he` is not supplied to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
the Hessian is obtained by differencing the gradient. With an analytic
gradient that is one numerical differentiation. When the gradient is
*also* differenced, the Hessian is a difference of differences, which
happens only in `newton()` and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md):
one Newton iteration then costs about \\4p^2\\ evaluations of the
objective, against about \\2p\\ evaluations of the gradient with an
analytic gradient. For an objective without either derivative,
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
which needs no Hessian, is usually cheaper.

## References

Dennis, J. E. and Schnabel, R. B. (1983). *Numerical Methods for
Unconstrained Optimization and Nonlinear Equations*. Prentice-Hall,
Englewood Cliffs, NJ.

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
# the start has eigenvalues 2 and -1.99; the repair is recorded in the
# trace.
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

# summary() counts the safeguards that fired during the run.
summary(minimize(newton(keep_trace = TRUE), rosen, c(-1.2, 1),
                 gr = rosen_gr, he = rosen_he))
#> <optimizer_result> Newton
#>   value      : 3.74398e-21
#>   par        : 1 1
#>   iterations : 21   evaluations: f 29, g 22
#>   elapsed    : 1 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
#>   safeguards :
#>     step shortened: 4
```
