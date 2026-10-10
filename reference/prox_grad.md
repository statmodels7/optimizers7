# Proximal Gradient Method

Minimizes \\f(x) + g(x)\\ where \\f\\ is smooth and \\g\\ is handled
entirely through its proximal operator, so that a non-differentiable
term is minimized without being differenced. With `accelerate = TRUE`
the method is an accelerated iteration of the kind of Beck and Teboulle,
with momentum weight \\(k-1)/(k+2)\\, whose objective gap falls as
\\O(1/k^2)\\ against the \\O(1/k)\\ of the plain iteration.

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
  allowance for its rounding error? Defaults to `TRUE`, and is ignored
  when `accelerate` is `FALSE`. See the details about a tight tolerance.

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
Box bounds are rejected: pass the constraint through `prox` instead,
composing the projection onto the box with the operator.

## Details

Each iteration takes a gradient step on the smooth part and applies the
proximal operator to the result, \$\$x\_{k+1} = \mathrm{prox}\_{t
g}\big(y_k - t\nabla f(y_k)\big),\$\$ with \\y_k = x_k\\ for the plain
method and \\y_k = x_k + \frac{k-1}{k+2}(x_k - x\_{k-1})\\ for the
accelerated one. The step length is found by backtracking: \\t\\ is
multiplied by `shrink` (one half by default) until the quadratic model
built at \\y_k\\ dominates \\f\\ at the new point, which is the
condition used by the convergence proof and needs no knowledge of the
Lipschitz constant. The step length is not increased again during the
run.

## What the stopping rule reads

The gradient of \\f + g\\ does not vanish at the solution, since \\g\\
is not differentiable there, so this method reports the **proximal
gradient mapping** \$\$G_t(x) = \frac{x - \mathrm{prox}\_{tg}(x -
t\nabla f(x))}{t}\$\$ as its gradient, read at the iterate \\x\\, the
point that is reported. It vanishes exactly at a stationary point of
\\f + g\\ and reduces to \\\nabla f\\ when \\g\\ is absent, so
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
keeps its meaning and its default tolerance. On the lasso of the
examples the mapping agrees with the KKT conditions, which the last
example checks.

Without momentum, the next iteration starts at the iterate, so the
gradient read there serves as the gradient of the next step, and the
plain method evaluates the gradient once per iteration. The accelerated
method takes its step at the extrapolated point, which is not the point
that it reports, and evaluates the gradient about twice per iteration.

## Acceleration

On an ill-conditioned problem the accelerated iteration needs far fewer
iterations than the plain one, as the \\O(1/k^2)\\ and \\O(1/k)\\ rates
indicate. On a well-conditioned problem with a tight tolerance the plain
iteration can need fewer: near the solution the objective changes at the
level of its rounding error, and each restart triggered there discards
the momentum built since the previous one. At the default
`crit_grad(1e-6)` this does not occur on the lasso of the examples.

## Restarting

Momentum makes the objective non-monotone, and an increase far from the
solution indicates momentum built in the wrong direction. With
`restart = TRUE` an increase resets the extrapolation to the current
point. This is the adaptive restart of O'Donoghue and Candes, and it
costs one comparison per iteration. On a badly conditioned problem it
reduces the number of iterations substantially.

An increase is measured against an allowance for the rounding error of
the objective, \\8\\\varepsilon\\\max(1, \lvert f \rvert)\\, where \\f\\
is the current total objective, and not against zero. The objective is a
sum, so its rounding error grows with the number of terms, and a bare
`>` test fires on that error once the iteration is close enough to the
solution, restarting the momentum at nearly every iteration. A larger
allowance would suppress real restarts on an ill-conditioned problem.
The error of a sum of \\n\\ terms can reach about \\n/2\\ units in the
last place, so for a large sum the allowance is smaller than the
possible error.

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
# iteration: stationarity where a coefficient is not zero, and a gradient
# inside the interval opened by the kink where it is.
gr_at <- gf(fit@par)
active <- fit@par != 0
max(abs(gr_at[active] + lambda * sign(fit@par[active])))
#> [1] 4.452503e-07
max(abs(gr_at[!active])) < lambda
#> [1] TRUE

# Box bounds are rejected, and the message names the argument that takes
# the constraint.
try(minimize(prox_grad(prox = function(v, t) v, g = function(b) 0),
             f, rep(0, 8), gr = gf, lower = 0))
#> Error : prox_grad() takes its constraint through the proximal operator, not through box bounds:
#>   compose the projection onto the box into 'prox'.
```
