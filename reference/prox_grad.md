# Proximal Gradient Method

Minimizes \\f(x) + g(x)\\ where \\f\\ is smooth and \\g\\ is handled
entirely through its proximal operator, so that a non-differentiable
term is minimized without ever differencing it. With `accelerate = TRUE`
the method is the accelerated one of Beck and Teboulle, whose objective
gap falls as \\O(1/k^2)\\ against the \\O(1/k)\\ of the plain iteration.

## Usage

``` r
prox_grad(
  prox,
  g,
  accelerate = TRUE,
  step = 1,
  shrink = 0.5,
  restart = TRUE,
  criterion = crit_grad(),
  maxit = 1000,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 20,
  keep_trace = FALSE
)
```

## Arguments

- prox:

  The proximal operator of the non-smooth part, a function of the point
  and the step length, `prox(v, step)`, returning the minimizer of
  \\\lVert b - v \rVert^2/(2\\\mathrm{step}) + g(b)\\.
  `penalties7::penalty_prox()` supplies one for every penalty that has
  it.

- g:

  The value of the non-smooth part, a function of the point. Required
  alongside `prox`: the two describe the same term, and without `g` the
  reported objective would be the smooth part alone.

- accelerate:

  Apply the momentum extrapolation? Defaults to `TRUE`.

- step:

  The initial step length offered to the backtracking search. Defaults
  to `1`.

- shrink:

  The factor by which a rejected step is reduced. Defaults to `0.5`.

- restart:

  Reset the momentum when the objective increases by more than the
  objective's own rounding? Defaults to `TRUE`, and is ignored when
  `accelerate` is `FALSE`. See the measured cost at a tight tolerance
  above.

- criterion:

  The stopping rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object. Defaults to
  [`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md),
  which here reads the proximal gradient mapping.

- maxit:

  Maximum iterations, a finite number at least 1. Defaults to 1000.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 20.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

## Value

An S7 object of class
[ProxGrad](https://statmodels7.github.io/optimizers7/reference/ProxGrad-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
Box bounds are refused: pass the constraint through `prox` instead,
composing the projection onto the box into it.

## Details

Each iteration takes a gradient step on the smooth part and applies the
proximal operator to the result, \$\$x\_{k+1} = \mathrm{prox}\_{t
g}\big(y_k - t\nabla f(y_k)\big),\$\$ with \\y_k = x_k\\ for the plain
method and \\y_k = x_k + \frac{k-1}{k+2}(x_k - x\_{k-1})\\ for the
accelerated one. The step length is found by backtracking: \\t\\ is
halved until the quadratic model built at \\y_k\\ dominates \\f\\ at the
new point, which is the condition the convergence proof uses and which
needs no knowledge of the Lipschitz constant.

## What the stopping rule reads

The gradient of \\f + g\\ does not vanish at the solution, \\g\\ being
non-differentiable there, so this method reports the **proximal gradient
mapping** \$\$G_t(x) = \frac{x - \mathrm{prox}\_{tg}(x - t\nabla
f(x))}{t}\$\$ as its gradient, read at the iterate. It vanishes exactly
at a stationary point of \\f + g\\ and reduces to \\\nabla f\\ when
\\g\\ is absent, so
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
keeps its meaning and its default tolerance. Measured on the lasso of
the examples below, the mapping agrees with the KKT violation: at the
reported point the active coordinates satisfy \\\lvert \nabla f +
\lambda \operatorname{sign}(\beta) \rvert \le\\ `4.5e-07` and the
inactive ones have \\\lvert \nabla f \rvert \le 0.129\\ against a
\\\lambda\\ of 0.4.

Acceleration pays 2.1 gradient evaluations per iteration against the
plain method's 1.0, the extrapolated point at which it takes its step
not being the point it reports.

## Acceleration, and where it earns its keep

On an ill-conditioned problem the difference is the one the theory
predicts. Measured on a smooth quadratic in eight unknowns at
`crit_grad(1e-8)`, plain against accelerated: 64 iterations against 31
at a condition number near 3, 1154 against 100 at 55, and 9321 against
268 at 480.

On a well-conditioned problem asked for a tight tolerance the plain
iteration is still the better one. The lasso below converges in 12
iterations at `crit_grad(1e-9)` with `accelerate = FALSE`, in 89 with
`restart = FALSE`, and in 782 with the defaults, all three reaching the
same support, the same objective to the last bit and the same
coefficients to `8e-09`. What costs the iterations is the restart, not
the momentum: near the solution the objective moves at the rounding
level, and each spurious reset discards the momentum built since the
last one. At the default `crit_grad(1e-6)` none of this appears, both
settings taking about ten iterations.

## Restarting

Momentum makes the objective non-monotone, and an increase far from the
solution is a symptom of momentum built in the wrong direction. With
`restart = TRUE` an increase resets the extrapolation to the current
point. This is the adaptive restart of O'Donoghue and Candes and costs
one comparison per iteration. It is worth a great deal where the problem
is badly conditioned: on a smooth quadratic in eight unknowns at
`crit_grad(1e-8)`, with the restart against without, 150 iterations
against 858 at a condition number of 55, 573 against 5418 at 480, and
1042 against 15104 at 2400.

An increase is measured against the objective's own rounding, not
against zero. The objective is a sum, so its error grows with the number
of terms, and a test reading a bare `>` fires on that error once the
iteration is near enough to the solution: the lasso below took 21646
iterations at `crit_grad(1e-9)` before the allowance and takes 782 after
it, at the same answer. The allowance is eight units in the last place
of the current objective, measured to give the same run anywhere between
one and thirty-two while 256 begins costing iterations on the
ill-conditioned quadratic, where it suppresses restarts that are real.
What it costs there is about a tenth: 149, 571 and 934 iterations at the
three condition numbers before, against the 150, 573 and 1042 above.

## References

Beck, A. and Teboulle, M. (2009). A fast iterative
shrinkage-thresholding algorithm for linear inverse problems. *SIAM
Journal on Imaging Sciences*, 2(1), 183–202.

O'Donoghue, B. and Candes, E. (2015). Adaptive restart for accelerated
gradient schemes. *Foundations of Computational Mathematics*, 15(3),
715–732.

## See also

[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
for a non-smooth method that needs no proximal operator,
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) for
the smooth case.

## Examples

``` r
# A lasso-penalized least squares problem, solved through the operator.
# Three of the eight coefficients are non-zero in the truth.
set.seed(1)
X <- matrix(rnorm(200 * 8), 200, 8)
b0 <- c(2, -1.5, 0, 0, 0.8, 0, 0, 0)
y <- as.numeric(X %*% b0 + rnorm(200))
lambda <- 0.4

f  <- function(b) sum((y - X %*% b)^2) / (2 * nrow(X))
gf <- function(b) as.numeric(-crossprod(X, y - X %*% b) / nrow(X))

fit <- minimize(
  prox_grad(prox = function(v, t) sign(v) * pmax(abs(v) - t * lambda, 0),
            g = function(b) lambda * sum(abs(b))),
  f, gr = gf, par = rep(0, 8))
round(fit@par, 3)
#> [1]  1.541 -1.091  0.000  0.000  0.420  0.000  0.000  0.000
sum(fit@par != 0)          # the operator sets coefficients exactly to zero
#> [1] 3

# The KKT conditions confirm the answer, and share no arithmetic with the
# iteration: stationary where a coefficient survives, and inside the
# interval the kink opens where one does not.
gr_at <- gf(fit@par)
active <- fit@par != 0
max(abs(gr_at[active] + lambda * sign(fit@par[active])))
#> [1] 4.452503e-07
max(abs(gr_at[!active])) < lambda
#> [1] TRUE

# Box bounds are refused, and the message says where the constraint goes.
try(minimize(prox_grad(prox = function(v, t) v, g = function(b) 0),
             f, rep(0, 8), gr = gf, lower = 0))
#> Error : prox_grad() takes its constraint through the proximal operator, not through box bounds:
#>   compose the projection onto the box into 'prox'.
```
