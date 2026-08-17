# Newton's Method with a Modified Hessian

Solves \\H d = -g\\ for the direction, repairing \\H\\ when it is not
positive definite, and then line searches along it.

## Usage

``` r
newton(
  criterion = crit_any(crit_grad(), crit_rel_obj()),
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

  The stopping rule; see
  [`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md).

- hessian_mod:

  How to repair an indefinite Hessian: `"eigen"` (default) or `"ridge"`.
  See Details.

- floor:

  The smallest eigenvalue the repaired Hessian is allowed. Defaults to
  `1e-8`.

- step:

  Initial step length offered to the line search. Defaults to `1`, which
  is the natural Newton step and is accepted unchanged near the
  solution.

- line_search:

  See
  [`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
  and
  [`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md).

- maxit, max_eval, verbose, refresh, keep_trace:

  As in
  [`gd`](https://statmodels7.github.io/optimizers7/reference/gd.md).

## Value

An S7 object of class `Newton`.

## Details

Newton's method converges quadratically near a minimum and is unreliable
away from one: \\H^{-1}g\\ is only a descent direction when \\H\\ is
positive definite. Where the objective curves downwards the unmodified
step points towards a saddle or a maximum, and no line search can rescue
it — every step along an ascent direction increases the objective. This
is not an edge case; it is the ordinary situation far from the solution.

Both repairs begin with an attempted Cholesky factorization, which when
it succeeds is simultaneously the test for positive definiteness and the
solve. When it fails:

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

Which repair fired is recorded in the trace, so the trace shows a run
that spent its time repairing rather than converging.

**On the Hessian itself.** If `he` is not supplied to
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
the Hessian is obtained by differencing the gradient. That is one
numerical differentiation when the gradient is analytic and acceptable;
when the gradient is *also* differenced it is two composed, which is the
one place in the package where that happens and where the result is
correspondingly poor. An objective with neither derivative is better
served by
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
which never needs a Hessian at all.

## References

Gill, P. E., Murray, W. and Wright, M. H. (1981). *Practical
Optimization*. Academic Press, London.

Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*, 2nd
edition. Springer, New York.

## See also

[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md),
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

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
#>   elapsed    : 1 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-12 (relative))
```
