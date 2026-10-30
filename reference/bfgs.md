# BFGS

Builds an approximation to the inverse Hessian from successive
gradients, so the direction is a matrix-vector product and the method
never evaluates second derivatives. It is the usual choice for a smooth
problem of moderate size.

## Usage

``` r
bfgs(
  criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
  curvature_tol = 1e-10,
  max_skip = 5,
  step = 1,
  line_search = wolfe(),
  maxit = 500,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 10,
  keep_trace = FALSE,
  stop_decrement = 0.001
)
```

## Arguments

- criterion:

  The stopping rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object. Defaults to
  `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`.

- curvature_tol:

  The update is skipped when \\s^\top y \le \texttt{curv\\tol}\\\lVert
  s\rVert\lVert y\rVert\\. Defaults to `1e-10`. A larger value skips
  more often. Zero is admitted and gives the textbook condition \\s^\top
  y \> 0\\. A negative value is rejected, because it would admit pairs
  with negative curvature, and the approximation would then no longer be
  guaranteed to stay positive definite.

- max_skip:

  Consecutive skipped updates before the approximation is reset to the
  identity. Defaults to 5, and has to be a whole number. Zero is
  admitted and resets on the first skip; a negative value is rejected,
  because the comparison `>=` would treat it as zero.

- step, line_search, maxit, max_eval, verbose, refresh, keep_trace:

  As in
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
  except that `line_search` defaults to
  [`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
  and `maxit` to 500. See below for the reason for the line search.

- stop_decrement:

  The largest decrease \\\tfrac12\lvert g^\top s\rvert\\, with \\s\\ the
  last step and \\g\\ the gradient at the new point, at which a stopping
  rule on the change in the parameters or in the objective, such as
  [`crit_abs_par()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md)
  or
  [`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md),
  may end the run as converged. A single positive number, in the units
  of the objective, or `Inf` to accept every such stop. Defaults to
  `1e-3`. See below.

## Value

An S7 object of class
[Bfgs](https://statmodels7.github.io/optimizers7/reference/Bfgs-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Notation

\\s = x\_{new} - x\_{old}\\ is the secant vector and \\y = g\_{new} -
g\_{old}\\ the change in the gradient. The scalar step length that a
line search returns is \\\alpha\\, a different quantity; the two are
related by \\s = \alpha d\\ for the direction \\d\\.

## Why the line search is the strong Wolfe one

The update is meaningful only when \\s^\top y \> 0\\, and the Wolfe
curvature condition guarantees that it holds: Armijo backtracking can
accept a step so short that the gradient has barely changed, which gives
a pair without curvature.

On the problems of
[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
the difference is small in practice. Under
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md),
BFGS skips an update on a few problems, needs a few more iterations on
two of them, and converges on the non-smooth `abs_sum`, where it stops
without converging under
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md).
Wolfe is the default because it guarantees a sound update; Armijo is a
reasonable choice where evaluations are expensive.

## Skipping and resetting

When the curvature condition fails, the update is skipped. A small
\\s^\top y\\ makes \\\rho = 1/s^\top y\\ very large, and one such update
can destroy the accumulated approximation; skipping it keeps the
previous matrix, which is still positive definite. After `max_skip`
consecutive skips the matrix is reset to the identity, so the method
restarts as steepest descent and rebuilds the approximation. Both events
appear in the trace, as `bfgs update skipped` and `bfgs reset`. After a
reset of this kind the next direction is \\-g\\, without the scaling
described below; after a reset caused by a direction that is not a
descent direction, the direction is scaled.

The first accepted pair rescales the identity by \\s^\top y / y^\top
y\\. Without it the first quasi-Newton step is taken with a unit
Hessian, which on a badly scaled problem has the wrong magnitude and
costs backtracks. Until a pair has been accepted, the direction is the
negative gradient multiplied by \\\min(1, 1/\lVert g\rVert\_\infty)\\,
which gives a displacement of order one in the parameters and never
lengthens the step.

## A stop on a change

A stopping rule on the change in the parameters or in the objective
states that the iterate no longer moves, and near a minimum that is the
same as being there. Away from one it is not: where the objective is
steep, or the curvature has been learned along a few directions only, a
step can be short while the objective is still falling. A stop on the
change is therefore not accepted where the decrease \\\tfrac12\lvert
g^\top s\rvert\\ that the last step still predicts at the new point
exceeds `stop_decrement`. The run ends as not converged, with
`criterion_met` reporting `"stopped without converging"` and the message
giving the decrease, and the trace records the safeguard
`change stop refused`. A rule on the gradient, such as
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md),
is not affected. The check uses quantities that the loop already holds.
`stop_decrement = Inf` accepts every stop on a change.

## References

Broyden, C. G. (1970). The convergence of a class of double-rank
minimization algorithms. *IMA Journal of Applied Mathematics* **6**,
76–90. The update was obtained independently the same year by Fletcher,
Goldfarb and Shanno, whence the name.

Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*, 2nd
edition. Springer, New York.

## See also

[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
for many parameters,
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
when a Hessian is available,
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
for the default line search of this method.

## Examples

``` r
rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
rg <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                    200 * (p[2] - p[1]^2))
minimize(bfgs(), rosen, c(-1.2, 1), gr = rg)
#> <optimizer_result> BFGS
#>   value      : 3.22826e-13
#>   par        : 1 1
#>   iterations : 35   evaluations: f 49, g 39
#>   elapsed    : 2 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)

# With Armijo instead of Wolfe, the update is skipped once here, the trace
# records it, and the run converges.
a <- minimize(bfgs(line_search = armijo(), keep_trace = TRUE), rosen,
              c(-1.2, 1), gr = rg)
c(iterations = a@iterations, skips = sum(a@trace$safeguard == "bfgs update skipped"))
#> iterations      skips 
#>         39          1 

# A curvature threshold that every pair fails to meet forces every update
# to be skipped; the run exhausts its budget, and the trace names both
# safeguards.
b <- minimize(bfgs(curvature_tol = 1e10, keep_trace = TRUE), rosen, c(-1.2, 1),
              gr = rg)
c(converged = b@converged, value = b@value)
#>    converged        value 
#> 0.000000e+00 4.354898e-06 
unique(b@trace$safeguard)
#> [1] "bfgs update skipped" "bfgs reset"         
```
