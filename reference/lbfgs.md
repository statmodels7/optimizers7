# Limited-Memory BFGS

Applies the BFGS update without forming the matrix: only the last
`memory` secant pairs are kept, and the direction comes from the
two-loop recursion. It costs \\O(mp)\\ in time and memory where
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
costs \\O(p^2)\\. The direction equals the one of the explicitly
assembled inverse built from the same pairs and the same initial
scaling. The two methods do not give the same directions in a run,
because
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
keeps every pair and scales its initial matrix by the first accepted
pair, whereas the recursion uses the last `memory` pairs and the most
recent one for the scaling.

## Usage

``` r
lbfgs(
  criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
  memory = 10,
  curvature_tol = 1e-10,
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

- memory:

  How many secant pairs to keep, a single positive whole number.
  Defaults to 10. A fractional value is rejected instead of truncated,
  since the count is passed to the compiled code as an integer.

- curvature_tol:

  A pair is not stored when \\s^\top y \le \texttt{curv\\tol}\\\lVert
  s\rVert\lVert y\rVert\\. Defaults to `1e-10`. Zero is admitted and
  gives \\s^\top y \> 0\\; a negative value is rejected, for the reason
  given on the page of
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md).

- step, line_search, maxit, max_eval, verbose, refresh, keep_trace:

  As in
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md).

- stop_decrement:

  The largest decrease \\\tfrac12\lvert g^\top s\rvert\\ at which a
  stopping rule on the change in the parameters or in the objective may
  end the run as converged, as in
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md).
  A single positive number, in the units of the objective, or `Inf` to
  accept every such stop. Defaults to `1e-3`. See below.

## Value

An S7 object of class
[Lbfgs](https://statmodels7.github.io/optimizers7/reference/Lbfgs-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Notation

\\s = x\_{new} - x\_{old}\\ is the secant vector and \\y = g\_{new} -
g\_{old}\\ the change in the gradient. The scalar step length is
\\\alpha\\, a different quantity.

## The saving

The two-loop recursion returns exactly the product of the explicitly
assembled inverse built from the same pairs and the same scaling (the
last example). The methods differ in cost: one iteration costs
\\O(p^2)\\ for
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and \\O(mp)\\ for `lbfgs()`, while the numbers of iterations of the two
methods on a dense quadratic are close. The difference in time becomes
large for a few hundred parameters, and for several thousand the
\\O(p^2)\\ storage and work of the full matrix become prohibitive.

The recursion is scaled at each iteration by \\s^\top y / y^\top y\\ of
the most recent pair. This scalar estimates the inverse curvature along
the latest step and takes the place of the initial matrix of the full
update.

## Choosing memory

A larger `memory` is not always better. Old pairs describe the curvature
at points that the iterate has left, and each pair costs \\O(p)\\ per
iteration. Beyond a handful of pairs the number of iterations changes
little while the arithmetic keeps rising; ten is the conventional
choice.

## A stop on a change

A stop on the change in the parameters or in the objective is not
accepted where the decrease \\\tfrac12\lvert g^\top s\rvert\\ that the
last step still predicts at the new point exceeds `stop_decrement`; the
run then ends as not converged, as described in
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md).
`stop_decrement = Inf` accepts every such stop.

## References

Nocedal, J. (1980). Updating quasi-Newton matrices with limited storage.
*Mathematics of Computation* **35**, 773–782.

Liu, D. C. and Nocedal, J. (1989). On the limited memory BFGS method for
large scale optimization. *Mathematical Programming* **45**, 503–528.

## See also

[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
for a few parameters,
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
when a Hessian is available.

## Examples

``` r
rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
rg <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                    200 * (p[2] - p[1]^2))
minimize(lbfgs(memory = 5), rosen, c(-1.2, 1), gr = rg)
#> <optimizer_result> L-BFGS
#>   value      : 1.94698e-14
#>   par        : 1 1
#>   iterations : 36   evaluations: f 49, g 40
#>   elapsed    : 3 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)

# The two-loop recursion is not an approximation to the full update: on
# the same pairs it returns the same product, to rounding error.
set.seed(1)
p <- 6
S <- matrix(rnorm(p * 4), p, 4)
Y <- matrix(rnorm(p * 4), p, 4) + 3 * S      # positive curvature
g <- rnorm(p)
gamma <- sum(S[, 4] * Y[, 4]) / sum(Y[, 4]^2)

H <- diag(gamma, p)                          # the explicit inverse
for (j in 1:4) {
  s <- S[, j]; y <- Y[, j]; rho <- 1 / sum(s * y)
  V <- diag(p) - rho * outer(s, y)
  H <- V %*% H %*% t(V) + rho * outer(s, s)
}

two_loop <- function(g, S, Y, gamma) {       # the recursion
  m <- ncol(S); a <- numeric(m); q <- g
  for (j in m:1) {
    rho <- 1 / sum(S[, j] * Y[, j])
    a[j] <- rho * sum(S[, j] * q); q <- q - a[j] * Y[, j]
  }
  r <- gamma * q
  for (j in 1:m) {
    rho <- 1 / sum(S[, j] * Y[, j])
    b <- rho * sum(Y[, j] * r); r <- r + S[, j] * (a[j] - b)
  }
  r
}
max(abs(H %*% g - two_loop(g, S, Y, gamma)))
#> [1] 1.665335e-16
```
