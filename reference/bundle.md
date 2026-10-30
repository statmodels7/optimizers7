# The Proximal Bundle Method

Proximal bundle method for convex non-smooth objectives: subgradients
collected at the visited points define a cutting-plane model, a proximal
term keeps the step near the stability center, and convergence is tested
on the optimality estimate \\\lVert p \rVert^2 + \alpha\\, built from
the aggregate subgradient \\p\\ and its linearization error \\\alpha\\,
and not on a single subgradient, which does not vanish at a kink.

## Usage

``` r
bundle(
  criterion = crit_stationary(),
  t0 = 1,
  t_min = 1e-10,
  t_max = 1e+10,
  accept_fraction = 0.1,
  bundle_size = 20,
  qp_iters = 500,
  qp_tol = 1e-12,
  maxit = 500,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 10,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule. Defaults to `crit_stationary(1e-8)`, read on the
  optimality estimate \\\lVert p \rVert^2 + \alpha\\.

- t0:

  Initial value of the proximal parameter \\t\\, expressed as a step
  length on the parameter scale and not as a bare multiplier; see
  Details. Defaults to `1`.

- t_min, t_max:

  Bounds on it, so that a run of null steps or a run of serious ones
  cannot drive it to zero or to infinity. The defaults are `1e-10` and
  `1e10`.

- accept_fraction:

  Fraction of the predicted decrease that a step must actually deliver
  to be accepted, in \\(0, 1)\\. Defaults to `0.1`.

- bundle_size:

  Largest number of linearizations kept; when a new one would exceed it,
  the oldest are replaced by the aggregate. Defaults to `20`.

- qp_iters, qp_tol:

  Effort spent on the subproblem: at most this many accelerated
  projected-gradient steps, stopping when the weights move by less than
  `qp_tol`. The defaults are `500` and `1e-12`.

- maxit:

  Maximum iterations. Defaults to 500.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`, which sets no
  evaluation budget, so the run stops on the criterion or on `maxit`.
  Set a finite value to cap the cost of a run.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 10.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

## Value

An S7 object of class
[Bundle](https://statmodels7.github.io/optimizers7/reference/Bundle-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Behavior at a kink

Near the minimum of \\\lvert x \rvert\\ the evaluated subgradient is
\\\pm 1\\. The gradient norm seen by a descent method therefore stays
large, the line search finds no acceptable step, and the run stops at
the minimizer without meeting its stopping rule, whatever the tolerance.
On the median example below,
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
reaches the median and reports `converged = FALSE`.

The bundle method tests no single subgradient. It keeps a *collection*
of them, from the points that it has visited, and builds the
piecewise-linear model \\\max_j \\ f_j + g_j'(x - x_j) \\\\, the largest
of the linearizations. Such a model represents a kink well, as the
meeting point of two linearizations. At the minimum of \\\lvert x
\rvert\\ the bundle holds subgradients near \\+1\\ and near \\-1\\, and
\\0\\ is a convex combination of them, which is the statement that \\0\\
lies in the subdifferential.

The step solves \$\$\min_d\\ \max_j \\ -\alpha_j + g_j'd \\ +
\frac{1}{2t}\lVert d \rVert^2,\$\$ where \\\alpha_j \ge 0\\ measures how
far linearization \\j\\ lies below the objective at the current point,
and the quadratic term keeps the step inside the region where the model
is trusted. The minimizer is \\d = -t\\p\\, with \\p\\ the aggregate
subgradient and \\\alpha\\ the aggregate linearization error. The value
of the piecewise-linear model at \\d\\, \\v = -(t\lVert p \rVert^2 +
\alpha) \le 0\\, is the *predicted decrease*, which the acceptance test
below reads.

[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
reads a related but different quantity, the *optimality estimate*
\\\lVert p \rVert^2 + \alpha\\. Both vanish at a solution, but \\-v\\
contains the factor \\t\\, and \\t\\ is halved at every null step, so
\\-v\\ can fall below any tolerance because the trust parameter has
shrunk and not because the point is stationary. A run that stopped on
\\-v\\ would then report convergence at a point where the model still
has a large aggregate subgradient. The optimality estimate omits the
factor \\t\\ and does not depend on the size of the trust parameter.

## Serious steps and null steps

A trial point is accepted, making a *serious step*, only if the decrease
that it delivers is at least the fraction `accept_fraction` of the
predicted decrease. This plays the role of the Armijo condition in a
line search: without a required fraction, a sequence of ever smaller
improvements could be accepted although it does not converge to a
minimizer.

A rejected trial point still contributes its subgradient to the bundle,
which improves the model at the next iteration. Such an iteration is a
*null step*, and a run with many null steps is refining the model of a
kink. Both counts appear in the `message` of the result, as in
`8 serious, 14 null`.

`t` is halved after a null step and doubled after a serious one, within
`t_min` and `t_max`. Kiwiel's rule chooses the factor from a curvature
estimate and performs better; the rule implemented here is simpler and
bounded.

## The trust parameter t0

The step is \\d = -t\\p\\, so with a bare `t` the length of the first
step would be proportional to the size of the gradient. On a steep
objective that step lands where the objective is enormous, each null
step halves `t` while the gradients collected there are larger still,
and the run can spend its whole budget on rejected steps.

When the largest absolute component of the gradient at the starting
point exceeds 1, `t0` is therefore divided by that component, so that
the largest coordinate of the first step has size `t0`; otherwise `t0`
is used as it is. This is the normalization of a line search that starts
at 1 along a unit direction, and `t_min` and `t_max` are scaled in the
same way. A larger `t0` suits a problem whose optimum is far away, and a
smaller one a problem on which the model is accurate only nearby.

Should a subgradient overflow anyway, which an objective that grows fast
enough can produce, the run stops with a message instead of spending its
budget on a subproblem whose matrix contains an infinity.

## Bounded memory

Without a limit the bundle grows indefinitely. When a new linearization
would make it exceed `bundle_size`, the oldest linearizations are
replaced by the *aggregate*, the affine function defined by the solution
of the subproblem at the current stability center. Discarding them would
lose their information and can stall the method; the aggregate retains a
summary of that information in one element.

## Convexity requirement

The theory is for convex \\f\\, where the linearization errors
\\\alpha_j\\ are non-negative automatically. On a non-convex objective
they can be negative and are clipped at zero. The method still runs and
usually works, but the convergence guarantee no longer applies: a
negative error means that \\f\\ lies below one of its own
linearizations, and setting the error to zero discards that information.

## Subgradients

Supply `gr`. Without it the package computes the gradient by finite
differences, which is accurate wherever \\f\\ is differentiable, that
is, almost everywhere. A difference taken *across* a kink is not a
subgradient, and the model is then built from a number that does not
correspond to any linearization of \\f\\.

## References

Kiwiel, K. C. (1990). Proximity control in bundle methods for convex
nondifferentiable minimization. *Mathematical Programming* **46**,
105–122.

Mäkelä, M. M. (2002). Survey of bundle methods for nonsmooth
optimization. *Optimization Methods and Software* **17**, 1–29.

## See also

[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
for the non-smooth case where a proximal operator is available,
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
and
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
for the derivative-free methods,
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for the default rule of this method.

## Examples

``` r
bundle()
#> <optimizer> proximal bundle
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 500, evaluations Inf
#>   settings  : t0 = 1, t_min = 1e-10, t_max = 1e+10, accept_fraction = 0.1, bundle_size = 20, qp_iters = 500, qp_tol = 1e-12

# The median is the minimizer of a sum of absolute deviations: the
# objective has a kink at every observation, and one of them is the answer.
set.seed(1)
y <- rnorm(101)
f <- function(p) sum(abs(y - p))
g <- function(p) -sum(sign(y - p))

r <- minimize(bundle(), f, par = 0, gr = g)
c(bundle = r@par, median = median(y))
#>     bundle     median 
#> 0.07456498 0.07456498 
abs(r@par - median(y))     # 8e-16
#> [1] 1.151856e-15
r@message                  # how many trials were accepted
#> [1] "8 serious, 14 null"

# The aggregate subgradient goes to zero at the kink, whereas the
# subgradient computed by g just beside it has norm 1.
r@gradient
#> [1] 2.767786e-13
abs(g(r@par - 1e-9))
#> [1] 1

# A descent method lands on the same point and cannot certify it.
d <- minimize(bfgs(criterion = crit_grad()), f, par = 0, gr = g)
c(par = d@par, converged = d@converged)
#>        par  converged 
#> 0.07456498 0.00000000 

# The kink lies along a diagonal. The bundle method and the simplex reach the
# minimum; the random-direction pattern search (the default of compass())
# comes close, whereas a coordinate search stops at a wrong point.
ff <- function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2)
gg <- function(p) c(sign(p[1] + p[2]), sign(p[1] + p[2])) + 0.2 * p
c(bundle      = minimize(bundle(), ff, c(1, 0.5), gr = gg)@value,
  nelder_mead = minimize(nelder_mead(), ff, c(1, 0.5))@value,
  compass     = minimize(compass(), ff, c(1, 0.5))@value)
#>       bundle  nelder_mead      compass 
#> 8.206869e-09 9.439361e-17 2.114785e-03 

# A gradient rule is rejected, since this method declares no gradient.
try(minimize(bundle(criterion = crit_grad()), f, par = 0, gr = g))
#> Error : The stopping rule needs gradient, which proximal bundle does not provide.
#>   Choose a criterion that this optimizer can evaluate, or a method that provides it.
```
