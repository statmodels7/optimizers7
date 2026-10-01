# Chebyshev's Third-Order Method

Newton's method with Chebyshev's correction: where the Hessian \\H\\ is
positive definite the Newton step \\d_N = -H^{-1}g\\ is corrected by the
third derivative of the objective, \$\$d = d_N - \tfrac12 H^{-1}\\
T\[d_N, d_N\],\$\$ where \\T\[u, v\]\_i = \sum\_{jk} \partial^3
f/\partial x_i\partial x_j \partial x_k\\ u_j v_k\\. Near a minimum the
iteration converges cubically where Newton's converges quadratically.
Everything else (the repair of an indefinite Hessian, the bound on the
length of a step, the line search, the stopping rule) is
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)'s.

## Usage

``` r
chebyshev(
  criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
  ratio = 0.5,
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

- criterion, hessian_mod, floor, step, line_search, maxit, max_eval,
  verbose, refresh, keep_trace, max_length, typical:

  As in
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

- ratio:

  The largest relative change of the Newton step the correction may
  make, \\\lVert d - d_N\rVert \le \texttt{ratio}\\\lVert d_N\rVert\\. A
  larger correction, or one that is not a descent direction, means the
  cubic model does not describe the objective over the length of the
  step, and the Newton step is used instead. Defaults to `0.5`.

## Value

An S7 object of class
[Chebyshev](https://statmodels7.github.io/optimizers7/reference/Chebyshev-class.md),
inheriting from
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)'s
class, to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## The third derivative

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
passes the method a function `t3(x, d)` returning the vector \\T\[d,
d\]\\, the third derivative contracted twice with \\d\\. The tensor
itself is never needed: for a sum over observations \\f = \sum_i
\phi_i(x_i^\top\beta)\\ it is \\T\[d, d\] = X^\top(\phi'''\\(Xd)^2)\\,
which costs as much as a gradient. Without `t3` the contraction is one
second difference of the gradient along \\d\\, at two gradient
evaluations an iteration.

## Where the correction is applied

Only at an iterate where \\H\\ is positive definite and its Newton step
is a descent direction; where the Hessian has to be repaired the step is
the repaired Newton step, as in
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).
The correction then costs one further solve with the same matrix. A
refused correction is recorded in the trace as
`cubic correction refused`.

Measured on the inner problems of statmodels7 (seven models, 14 to 113
coefficients) started from the mode at neighbouring hyperparameters, the
method reaches the mode in one iteration fewer than
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
on six of seven, at the same point. Far from the mode the correction is
refused in most iterations, and the method behaves as
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

Inside whole statmodels7 fits the gain disappears. Over the 32 models of
its reference battery on which every inner method reaches the same
point, the median number of inner iterations is 1.00 times
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)'s
and the median time 1.01 times: an inner fit restarted from the previous
mode is usually at its mode after one or two Newton steps, so a faster
local rate has nothing left to shorten. It reduces the iterations where
an inner fit starts far from its mode (761 against 534 on a lasso beside
a random effect) and is no use as a replacement for that package's own
scoring iteration.

## Why not Halley

Halley's method solves with \\H + \tfrac12 T\[d_N\]\\ in place of \\H\\,
which needs the contraction as a matrix and a second factorization, of a
matrix that need not be positive definite. On the same models it took
the same number of iterations as Chebyshev's and more time.

## Box constraints

With `lower` or `upper` the method runs on the unconstrained scale, as
every method here does. A `t3` written for the original parameters is
then the wrong derivative, so it is rejected together with bounds; the
second difference of the gradient serves there.

## References

Gundersen, G. and Steihaug, T. (2010). On large-scale unconstrained
optimization problems and higher order methods. *Optimization Methods
and Software* **25**, 337–358.

Gutierrez, J. M. and Hernandez, M. A. (1997). A family of
Chebyshev-Halley type methods in Banach spaces. *Bulletin of the
Australian Mathematical Society* **55**, 113–130.

## See also

[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
whose settings it shares.

## Examples

``` r
rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
rosen_gr <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                          200 * (p[2] - p[1]^2))
rosen_he <- function(p) matrix(
  c(2 - 400 * (p[2] - 3 * p[1]^2), -400 * p[1],
    -400 * p[1], 200), 2, 2)
# T[d, d] of Rosenbrock: the only third derivatives are
# f_111 = 2400 p1 and f_112 = -400
rosen_t3 <- function(p, d) c(2400 * p[1] * d[1]^2 - 800 * d[1] * d[2],
                             -400 * d[1]^2)

a <- minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
b <- minimize(chebyshev(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he,
              t3 = rosen_t3)
c(newton = a@iterations, chebyshev = b@iterations)
#>    newton chebyshev 
#>        21        18 

# Without t3 the contraction is a second difference of the gradient.
minimize(chebyshev(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
#> <optimizer_result> Chebyshev
#>   value      : 0
#>   par        : 1 1
#>   iterations : 18   evaluations: f 30, g 53
#>   elapsed    : 1 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
```
