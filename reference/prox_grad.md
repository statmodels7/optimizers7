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
  \\\lVert b - v \rVert^2/(2\\\mathrm{step}) + g(b)\\. `penalty_prox`
  supplies one for every penalty that has it.

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

  Reset the momentum when the objective increases? Defaults to `TRUE`,
  and is ignored when `accelerate` is `FALSE`.

- criterion:

  The stopping rule; see
  [`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md).

- maxit:

  Maximum iterations. Defaults to 1000.

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
[`ProxGrad`](https://statmodels7.github.io/optimizers7/reference/ProxGrad-class.md),
to be handed to
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Details

Each iteration takes a gradient step on the smooth part and applies the
proximal operator to the result, \$\$x\_{k+1} = \mathrm{prox}\_{t
g}\big(y_k - t\nabla f(y_k)\big),\$\$ with \\y_k = x_k\\ for the plain
method and \\y_k = x_k + \frac{k-1}{k+2}(x_k - x\_{k-1})\\ for the
accelerated one. The step length is found by backtracking: \\t\\ is
halved until the quadratic model built at \\y_k\\ dominates \\f\\ at the
new point, which is the condition the convergence proof uses and which
needs no knowledge of the Lipschitz constant.

### What the stopping rule reads

The gradient of the total objective does not vanish at the solution –
that is what non-differentiable means – so this method reports the
**proximal gradient mapping** \$\$G_t(x) = \frac{x -
\mathrm{prox}\_{tg}(x - t\nabla f(x))}{t}\$\$ as its gradient, read at
the iterate. It vanishes exactly at a stationary point of \\f + g\\ and
reduces to \\\nabla f\\ when \\g\\ is absent, so
[`crit_grad`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
keeps its meaning and its default tolerance. The accelerated variant
pays one extra gradient per iteration for it, the extrapolated point at
which it takes its step not being the point it reports.

### Restarting

Momentum makes the objective non-monotone, and an increase far from the
solution is a symptom of momentum built in the wrong direction. With
`restart = TRUE` an increase resets the extrapolation to the current
point, which is the adaptive restart of O'Donoghue and Candes and costs
one comparison per iteration.

## References

Beck, A. and Teboulle, M. (2009). A fast iterative
shrinkage-thresholding algorithm for linear inverse problems. *SIAM
Journal on Imaging Sciences*, 2(1), 183–202.

O'Donoghue, B. and Candes, E. (2015). Adaptive restart for accelerated
gradient schemes. *Foundations of Computational Mathematics*, 15(3),
715–732.

## See also

[`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
for a non-smooth method that needs no proximal operator,
[`gd`](https://statmodels7.github.io/optimizers7/reference/gd.md) for
the smooth case.

## Examples

``` r
# a lasso-penalized least squares problem, solved through the operator
set.seed(1)
X <- matrix(rnorm(200 * 8), 200, 8)
b0 <- c(2, -1.5, 0, 0, 0.8, 0, 0, 0)
y <- as.numeric(X %*% b0 + rnorm(200))
lambda <- 0.4

fit <- minimize(
  prox_grad(prox = function(v, t) sign(v) * pmax(abs(v) - t * lambda, 0),
            g = function(b) lambda * sum(abs(b))),
  fn = function(b) sum((y - X %*% b)^2) / (2 * nrow(X)),
  gr = function(b) -crossprod(X, y - X %*% b) / nrow(X),
  par = rep(0, 8))
round(fit@par, 3)
#> [1]  1.541 -1.091  0.000  0.000  0.420  0.000  0.000  0.000
```
