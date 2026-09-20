# BFGS

Builds an approximation to the inverse Hessian from successive
gradients, so the direction is a matrix-vector product and no second
derivatives are ever required. The default method for a smooth problem
of moderate size: on Rosenbrock from the customary start it converges in
35 iterations and 49 objective evaluations with the gradient supplied.

## Usage

``` r
bfgs(
  criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
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

  The stopping rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object. Defaults to
  `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`.

- curv_tol:

  The update is skipped when \\s^\top y \le \texttt{curv\\tol}\\\lVert
  s\rVert\lVert y\rVert\\. Defaults to `1e-10`. A larger value skips
  more often. Zero is admitted and gives the textbook condition \\s^\top
  y \> 0\\; a negative value is refused, because the comparison would
  then hold for every pair and the skip protection would be off without
  saying so.

- max_skip:

  Consecutive skipped updates before the approximation is reset to the
  identity. Defaults to 5, and has to be a whole number. Zero is
  admitted and resets on the first skip; a negative value is refused,
  the comparison being `>=`, which would make it behave as zero.

- step, line_search, maxit, max_eval, verbose, refresh, keep_trace:

  As in
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
  except that `line_search` defaults to
  [`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
  and `maxit` to 500. See below for why the line search matters here.

## Value

An S7 object of class
[Bfgs](https://statmodels7.github.io/optimizers7/reference/Bfgs-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Notation

\\s = x\_{new} - x\_{old}\\ is the **secant vector** and \\y =
g\_{new} - g\_{old}\\ the change in the gradient. The scalar step length
a line search chooses is \\\alpha\\, a different quantity; the two are
related by \\s = \alpha d\\ for the direction \\d\\.

## Why the line search is the strong Wolfe one

The update is meaningful only when \\s^\top y \> 0\\, and the Wolfe
curvature condition is exactly the guarantee that it holds: Armijo
backtracking can accept a step so short that the gradient has barely
moved, leaving a pair with no curvature in it.

In practice the guarantee is worth less than it sounds, and the honest
measurement is worth having. Over the eight problems of
[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md),
BFGS under Armijo skips an update on three of them, once each, and
matches Wolfe elsewhere: Rosenbrock 39 iterations against 35, himmelblau
11 against 10, and on the non-smooth `abs_sum` Armijo is the **better**
of the two, reaching `6.8e-08` in 37 iterations where Wolfe stops at
`1.3e-02` reporting failure. Wolfe is the default because it makes the
update sound by construction; Armijo is a reasonable choice where
evaluations are dear.

## Skipping and resetting

When the curvature condition fails anyway the update is **skipped**. A
small \\s^\top y\\ makes \\\rho = 1/s^\top y\\ enormous and one bad step
destroys the accumulated approximation; a stale but sound matrix beats a
fresh but corrupted one. After `max_skip` consecutive skips there is no
curvature information left worth keeping and the matrix is reset to the
identity, so the method restarts as steepest descent and rebuilds. Both
events appear in the trace, as `bfgs update skipped` and `bfgs reset`.

The first accepted pair rescales the identity by \\s^\top y / y^\top
y\\. Without it the first quasi-Newton step is taken with a unit
Hessian, which on a badly scaled problem has the wrong magnitude
entirely and wastes a line search discovering so. The first direction is
also scaled to \\\min(1, 1/\lVert g\rVert\_\infty)\\, a displacement of
order one in the parameters, which can only ever shorten it.

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
for the line search this method wants.

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
#>   elapsed    : 3 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)

# Armijo instead of Wolfe: the update is skipped once here, which the trace
# reports, and the run still arrives.
a <- minimize(bfgs(line_search = armijo(), keep_trace = TRUE), rosen,
              c(-1.2, 1), gr = rg)
c(iterations = a@iterations, skips = sum(a@trace$safeguard == "bfgs update skipped"))
#> iterations      skips 
#>         39          1 

# Forcing the skip: a curvature threshold nothing can meet exhausts the
# budget, and the trace names both safeguards.
b <- minimize(bfgs(curv_tol = 1e10, keep_trace = TRUE), rosen, c(-1.2, 1),
              gr = rg)
c(converged = b@converged, value = b@value)
#>    converged        value 
#> 0.000000e+00 4.354898e-06 
unique(b@trace$safeguard)
#> [1] "bfgs update skipped" "bfgs reset"         
```
