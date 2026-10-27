# Chebyshev's Third-Order Method

Newton's method with Chebyshev's correction: where the Hessian \\H\\ is
positive definite the Newton step \\d_N = -H^{-1}g\\ is corrected by the
third derivative of the objective, \$\$d = d_N - \tfrac12 H^{-1}\\
T\[d_N, d_N\],\$\$ where \\T\[u, v\]\_i = \sum\_{jk} \partial^3
f/\partial x_i\partial x_j \partial x_k\\ u_j v_k\\. Near a minimum the
iteration converges cubically where Newton's converges quadratically.
Everything else (the repair of an indefinite Hessian, the bound on the
length of a step, the line search, the stopping rule) is as in
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

## Usage

``` r
chebyshev(
  criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
  ratio = 0.5,
  modify = modify_eigen(),
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

- criterion, modify, step, line_search, maxit, max_eval, verbose,
  refresh, keep_trace, max_length, typical:

  As in
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).
  `modify` takes a
  [`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md)
  object
  ([`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md)
  by default, or
  [`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md));
  see
  [modifications](https://statmodels7.github.io/optimizers7/reference/modifications.md).

- ratio:

  The largest relative change of the Newton step that the correction may
  make, \\\lVert d - d_N\rVert \le \texttt{ratio}\\\lVert d_N\rVert\\. A
  larger correction, or one that is not a descent direction, indicates
  that the cubic model does not describe the objective over the length
  of the step, and the Newton step is used instead. Defaults to `0.5`.

## Value

An S7 object of class
[Chebyshev](https://statmodels7.github.io/optimizers7/reference/Chebyshev-class.md),
inheriting from
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)'s
class, to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## The third derivative

The caller supplies to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
a function `t3(x, d)` that returns the vector \\T\[d, d\]\\, the third
derivative contracted twice with \\d\\, and
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
passes it to the method. The tensor itself is never needed: for a sum
over observations \\f = \sum_i \phi_i(x_i^\top\beta)\\ the contraction
is \\T\[d, d\] = X^\top(\phi'''\\(Xd)^2)\\, which costs as much as a
gradient. Without `t3` the contraction is one second difference of the
gradient along \\d\\, at two gradient evaluations in each iteration at
which the correction is attempted.

## Where the correction is applied

The correction is attempted only at an iterate where \\H\\ is positive
definite and its Newton step is a descent direction; where the Hessian
has to be repaired, the step is the repaired Newton step, as in
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).
The correction costs one further solve with the same matrix. A
correction that is not used is recorded in the trace as
`cubic correction refused`.

Near the mode the correction is usually accepted, and the method needs
fewer iterations than
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).
Far from the mode the correction is usually rejected, and the method
then behaves as
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).
When the starting point is already within one or two Newton steps of the
mode, as for an inner fit restarted from the previous mode, a faster
local rate has little to shorten, and the method gives no gain over
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

## Halley's method

Halley's method solves with \\H + \tfrac12 T\[d_N\]\\ in place of \\H\\.
This needs the contraction as a matrix and a second factorization, of a
matrix that need not be positive definite, and the package does not
provide it.

## Box constraints

With `lower` or `upper` the method runs on the unconstrained scale, as
every method that accepts bounds does. A `t3` written for the original
parameters is then the wrong derivative, so `t3` together with bounds is
rejected; without `t3` the second difference of the gradient is used.

## References

Gundersen, G. and Steihaug, T. (2010). On large-scale unconstrained
optimization problems and higher order methods. *Optimization Methods
and Software* **25**, 337–358.

Gutierrez, J. M. and Hernandez, M. A. (1997). A family of
Chebyshev-Halley type methods in Banach spaces. *Bulletin of the
Australian Mathematical Society* **55**, 113–130.

## See also

[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
whose settings this method shares.

## Examples

``` r
rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
rosen_gr <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                          200 * (p[2] - p[1]^2))
rosen_he <- function(p) matrix(
  c(2 - 400 * (p[2] - 3 * p[1]^2), -400 * p[1],
    -400 * p[1], 200), 2, 2)
# The example computes T[d, d] of Rosenbrock: the only third derivatives are
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
#>   elapsed    : 5 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
```
