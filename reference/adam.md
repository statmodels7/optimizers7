# Adaptive Moment Estimation

Adaptive moment estimation (Adam) is a first-order method whose
coordinate-wise step lengths come from exponentially weighted first and
second moments of the gradient. It takes no line search and tolerates a
gradient that is correct only on average, which makes it the method
suited to a stochastic objective.

## Usage

``` r
adam(
  criterion = crit_never(),
  alpha = 0.01,
  beta1 = 0.9,
  beta2 = 0.999,
  eps = 1e-08,
  decay = 0,
  amsgrad = FALSE,
  maxit = 1000,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 100,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule. Defaults to
  [`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md),
  so the run is governed by `maxit`; see Details.

- alpha:

  The learning rate, a single positive number: the size of a step when
  the gradient is steady. Defaults to `0.01`.

- beta1:

  Decay rate of the first moment, the smoothed gradient. Defaults to
  `0.9`.

- beta2:

  Decay rate of the second moment, the smoothed squared gradient.
  Defaults to `0.999`.

- eps:

  Added to the square root of the second moment before dividing, so that
  a coordinate whose gradient has been uniformly zero does not divide by
  zero. Defaults to `1e-8`.

- decay:

  Reduces the learning rate as \\\alpha_t = \alpha/(1 + \delta (t -
  1))\\ at iteration \\t = 1, 2, \ldots\\, so that the first step uses
  \\\alpha\\. Defaults to `0`, a constant rate; see Details.

- amsgrad:

  Hold the second moment at its running maximum? Defaults to `FALSE`;
  see Details.

- maxit:

  Maximum iterations. Defaults to 1000, more than the Newton and
  quasi-Newton methods, because Adam takes many small steps where a
  second-order method takes few large ones.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`, which sets no
  evaluation budget, so the run stops on the criterion or on `maxit`.
  Set a finite value to cap the cost of a run.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 100.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

## Value

An S7 object of class
[Adam](https://statmodels7.github.io/optimizers7/reference/Adam-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Details

Adam keeps an exponentially weighted average of the gradient, \\m_t\\,
and of its square, \\v_t\\, and steps \$\$x\_{t+1} = x_t - \alpha\\ \hat
m_t / (\sqrt{\hat v_t} + \epsilon).\$\$ The division is elementwise: a
coordinate whose gradient has been consistently large is divided by a
large number and moves modestly, while one whose gradient is small but
persistent still moves. This is a diagonal preconditioner assembled from
the gradients already computed, so it costs nothing beyond them. A
diagonal cannot represent the correlation between parameters that
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
learns from the same information.

Both averages start at zero, so in the first iterations they are biased
towards zero; dividing by \\1 - \beta^t\\ removes that bias, without
which the first iterations would barely move.

## Relation to the descent methods

Adam takes no line search and does not require the objective to decrease
at every step. Because of this, Adam tolerates a gradient that is
correct only on average, and it provides no guarantee of monotone
progress: the run may end at a point worse than a point through which it
passed.

Adam is therefore a poor choice for a small smooth problem where a
Hessian is affordable: on a quadratic
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
needs a few iterations and Adam some hundreds (see the examples). Adam
suits a long parameter vector, a noisy objective, or a surface on which
a quadratic model is a poor description.

## Stochastic objectives

Adam draws no subsamples of its own. The optimizer does not represent
observations, and a stochastic objective is an ordinary function of the
parameters, so a minibatch is drawn inside the objective, as in the
following closure:

    batch <- function(par) {
      i <- sample.int(n, size = 0.05 * n)
      sum((y[i] - par)^2) / 2
    }
    minimize(adam(), batch, par = 0, gr = batch_gr)

Adam then behaves as on a minibatch of its own drawing, and
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs the run
because the draws happen in the code of the caller.

Resample inside the objective, not around the run. Calling
`minimize(adam(maxit = 1), ...)` in a loop with a new batch each time
does not work: \\m\\ and \\v\\ start at zero and the bias correction
restarts at \\t = 1\\, so every call takes a first step of length
\\\alpha\\, and the accumulated moments, on which the method depends,
are discarded at each call.

The gradient check that
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
makes before starting compares `fn` and `gr` at one point, and on a
resampling objective the two use different minibatches, so the check
warns on every stochastic fit, and the warning does not indicate an
error. Set `options(optimizers7.check_gradient = FALSE)` around such a
fit.

## Stopping

The default criterion is
[`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md),
so the run ends when `maxit` is reached and reports `converged = FALSE`,
because the run does not test any stopping rule. On a noisy objective
Adam with a fixed `alpha` circles an optimum instead of settling on it,
so a tolerance on the gradient seldom fires there.

On an exact objective a stopping rule such as
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
can be passed and fires. On a noisy objective the value and the gradient
are both estimates, so a rule that reads either is unreliable, and the
method does not distinguish an exact objective from a noisy one. The
rule is also tested at the point produced by the last iteration of the
budget.

## The safeguards

`eps` floors the denominator, so a coordinate whose gradient has been
uniformly zero is not divided by zero.

`decay` makes the learning rate \\O(1/t)\\, a schedule that satisfies
the Robbins-Monro conditions (a divergent sum and a convergent sum of
squares), which a run on a noisy objective needs in order to settle at
the optimum instead of fluctuating about it within a radius set by
\\\alpha\\. The minibatch example below compares the scatter of five
runs with and without decay. It is off by default, because an exact
objective has no noise to average away.

`amsgrad` replaces \\v_t\\ by its running maximum. Reddi, Kale and Kumar
(2018) exhibited a convex problem on which Adam as published fails to
converge, because \\v_t\\ can shrink and let a single large gradient
dominate the iterate long after it has passed. The maximum prevents
this, at the price that the second moment in the effective rate
\\\alpha/\sqrt{\max_s v_s}\\ of a coordinate can no longer decrease. It
is `FALSE` by default so that `adam()` is Adam as published; it can be
set to `TRUE` when a run fails to settle.

A gradient or an update that is not finite ends the run instead of
propagating a `NaN` into every later iterate. The result reports
`converged = FALSE`, with the message `gradient not finite; stopped` for
a gradient and `the update is not finite; stopped` for an update, and
`par` is the last usable point.

## References

Kingma, D. P. and Ba, J. (2015). Adam: A Method for Stochastic
Optimization. *ICLR*.

Reddi, S. J., Kale, S. and Kumar, S. (2018). On the Convergence of Adam
and Beyond. *ICLR*.

## See also

[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
for the smooth case,
[`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md)
for the default stopping rule,
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) for
the other method that accepts uphill steps by design.

## Examples

``` r
adam()
#> <optimizer> adam
#>   stop when : iteration budget
#>   budgets   : maxit 1000, evaluations Inf
#>   settings  : alpha = 0.01, beta1 = 0.9, beta2 = 0.999, eps = 1e-08, decay = 0, amsgrad = FALSE
adam(alpha = 0.05, amsgrad = TRUE)
#> <optimizer> adam
#>   stop when : iteration budget
#>   budgets   : maxit 1000, evaluations Inf
#>   settings  : alpha = 0.05, beta1 = 0.9, beta2 = 0.999, eps = 1e-08, decay = 0, amsgrad = TRUE

# On an exact quadratic it reaches the minimum, at the cost of many more
# evaluations than BFGS. The flag is FALSE because crit_never() tests
# nothing.
q <- function(p) sum((p - c(1, 2))^2)
qg <- function(p) 2 * (p - c(1, 2))
a <- minimize(adam(alpha = 0.1, maxit = 2000), q, c(0, 0), gr = qg)
c(a@par, evaluations = a@counts[["f"]], converged = a@converged)
#>                         evaluations   converged 
#>           1           2        2002           0 
minimize(bfgs(), q, c(0, 0), gr = qg)@counts[["f"]]
#> [1] 3

# A stochastic objective: the minibatch is drawn inside the function, so
# the optimizer needs no knowledge of it. The check that compares fn with
# gr would compare two different batches, so it is turned off.
set.seed(1)
y <- rnorm(2000, mean = 3)
m <- 100
batch    <- function(p) { i <- sample.int(2000, m); sum((y[i] - p)^2) / 2 }
batch_gr <- function(p) { i <- sample.int(2000, m); -sum(y[i] - p) }

old <- options(optimizers7.check_gradient = FALSE)
minimize(adam(alpha = 0.05, decay = 0.01, maxit = 2000),
         batch, par = 0, gr = batch_gr)@par
#> [1] 2.980926
mean(y)
#> [1] 2.986045

# A decaying rate settles the run: the scatter over five runs is smaller.
spread <- function(d) sd(replicate(5,
  minimize(adam(alpha = 0.05, decay = d, maxit = 2000),
           batch, par = 0, gr = batch_gr)@par))
set.seed(2); c(none = spread(0), some = spread(0.01))
#>        none        some 
#> 0.032751298 0.007941274 
options(old)

# A gradient that is not finite ends the run, and the message records it.
bad <- minimize(adam(maxit = 50), q, c(0, 0), gr = function(p) c(NaN, 1))
c(bad@converged, bad@message)
#> [1] "FALSE"                        "gradient not finite; stopped"
```
