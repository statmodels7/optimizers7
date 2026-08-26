# Limited-Memory BFGS

BFGS without ever forming the matrix: only the last `memory` secant
pairs are kept, and the direction comes from the two-loop recursion.
Costs \\O(mp)\\ in time and memory where
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
costs \\O(p^2)\\, and returns the same direction the full matrix would
when the pairs are the same.

## Usage

``` r
lbfgs(
  criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
  memory = 10,
  curv_tol = 1e-10,
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

- memory:

  How many secant pairs to keep, a single positive whole number.
  Defaults to 10. A fractional value is refused rather than truncated,
  the kernel reading the count through
  [`as.integer()`](https://rdrr.io/r/base/integer.html).

- curv_tol:

  A pair is discarded when \\s^\top y \le \texttt{curv\\tol}\\\lVert
  s\rVert\lVert y\rVert\\. Defaults to `1e-10`. Zero is admitted and
  gives \\s^\top y \> 0\\; a negative value is refused, for the reason
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)'s
  page gives.

- step, line_search, maxit, max_eval, verbose, refresh, keep_trace:

  As in
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md).

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

## Where the saving is

The two-loop recursion returns exactly the product the explicitly
assembled inverse would: built from the same four pairs and the same
scaling, the two agree to `1.7e-16`. What differs is the cost.

Measured on a dense quadratic, both methods taking the gradient and
converging to the same point:

|       |                |                |
|-------|----------------|----------------|
| **p** | **bfgs**       | **lbfgs**      |
| 5     | 9 iterations   | 8 iterations   |
| 50    | 15             | 15             |
| 200   | 16, 0.06 s     | 17, 0.00 s     |
| 800   | 16, **4.06 s** | 17, **0.02 s** |

The iteration counts barely differ; the cost per iteration is what
diverges, and it does so between 200 and 800 parameters. Beyond a few
thousand the full matrix is not storable at all.

The recursion is scaled at each iteration by the most recent pair's
\\s^\top y / y^\top y\\. That single number does the work the full
matrix would otherwise do, and is the reason the method converges
without one.

## Choosing memory

More is not better. Old pairs describe curvature at points the iterate
has left, and they cost \\O(p)\\ each per iteration. Measured on
Rosenbrock, `memory` of 3, 5, 10, 30 and 100 gives 38, 36, 37, 37 and 37
iterations: past a handful there is nothing left to gain and the
arithmetic keeps rising. Ten is the conventional choice for that reason.

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

# The two-loop recursion is not an approximation to the full update: on the
# same pairs it returns the same product, to machine precision.
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
