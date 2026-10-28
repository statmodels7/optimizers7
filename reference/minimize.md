# Minimize a Function

Runs an optimizer on an objective and returns the point that it reached,
the value there, and the rule that stopped it. This is the single entry
point: every algorithm in the package is a method of this generic, so a
change of method is a change of the first argument. Box constraints are
accepted by every method except
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md),
which takes its constraint through the proximal operator, and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
when it is given `t3`. Every method minimizes;
[`maximize()`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
negates the objective and restores the sign of the answer.

## Usage

``` r
minimize(
  optimizer,
  fn,
  par,
  gr = NULL,
  hess = NULL,
  lower = -Inf,
  upper = Inf,
  ...
)
```

## Arguments

- optimizer:

  An
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  object carrying the algorithm and its settings. Dispatch reads this
  argument only.

- fn:

  The objective: a function of the parameter vector returning a single
  number, to be minimized. An object of another class works whenever
  [`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md)
  has a method for it.

- par:

  A numeric vector of starting values, or a starter object; see Starters
  below. With bounds it must lie strictly inside them.

- gr:

  The gradient, a function of the parameter vector. `NULL`, the default,
  has it differenced from `fn`. It is ignored when `fn` is an objective
  that carries its own gradient.

- hess:

  The Hessian, a function of the parameter vector. `NULL` is the
  default.
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
  and
  [`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
  read it; every other method accepts it and does not use it, so calling
  code need not branch on the algorithm.

- lower, upper:

  Box constraints, numeric of length 1 (applying to every parameter) or
  of length `p` (one each). Any other length is an error. The defaults
  `-Inf` and `Inf` are no constraint, and a coordinate whose pair is
  infinite on both sides is left untransformed. `lower` must be strictly
  below `upper` in every coordinate.

- ...:

  Passed to the method selected by dispatch. Among the shipped methods,
  only the method for
  [`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
  reads from it, its argument `t3`.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
object: `par`, `value`, `gradient`, `counts`, `iterations`, `converged`,
`criterion_met`, `message`, `trace`, `optimizer`, `elapsed` and `seed`.

## Details

The problem solved is

\$\$\min\_{x \in \mathbb{R}^{p}} f(x) \qquad \text{subject to} \quad l
\le x \le u,\$\$

with \\f\\ the objective, \\l\\ and \\u\\ the bounds, and the
inequalities read coordinatewise. Every method reports the point where
it stopped together with the rule that stopped it. `converged` is `TRUE`
when the stopping rule is satisfied; it is also `TRUE` when a line
search with a declared resolution establishes that every step fails to
decrease the objective by more than that resolution, and when
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
reduces its bracket to the resolution of the iterate, and the result
then reports that stop with its own label in `criterion_met`. It is
never inferred from the end of the run, so a run that exhausts its
budget comes back with `converged = FALSE`.

Dispatch is on `optimizer` alone, so each algorithm is written once. The
objective is normalized separately, by
[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md),
which dispatches on `fn`; a caller with its own kind of objective
registers one method there and every algorithm accepts it.

## Derivatives, supplied and differenced

A gradient that is not supplied is computed by central finite
differences, at a cost of \\2p\\ evaluations of the objective each time,
and the `message` of the result records that the gradient was obtained
by finite differences.

A gradient that *is* supplied is checked before the run, with one
central difference along the gradient direction at `par`, which costs
two evaluations. If the difference disagrees with the gradient, it is
repeated at steps ten and a hundred times smaller, and a disagreement
that persists gives a warning that names both rates, as in

    'gr' does not appear to be the gradient of 'fn': along the gradient
    direction at 'par', 'fn' changes at rate 4.47 where 'gr' predicts 22.4.

and the run proceeds. Without the check, a `gr` computed from a
different model than `fn` appears only as a line-search failure at the
first iteration, with no indication of its cause. Set
`options(optimizers7.check_gradient = FALSE)` to turn the check off.

## Starters

`par` may be a starter object instead of a vector:
[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
for all zeros,
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md)
for a uniform draw from a chosen range. Both work on the unconstrained
scale and are mapped back through the bounds, so one constant gives a
sensible value for every kind of parameter (zero becomes one for a
variance, one half for a probability) and every draw lies inside its
box.

The number of parameters is obtained in one of three ways. It can be
given with `start_zeros(npar = 3)`. Otherwise a `lower` or `upper` of
more than one element gives it, since bounds are one per parameter.
Otherwise the objective is probed by
[`infer_npar()`](https://statmodels7.github.io/optimizers7/reference/infer_npar.md),
which tries the lengths 1 to 50 once, before the run begins, and
succeeds when exactly one of them is accepted.

The probe works for an objective with a fixed width built into it, such
as `X %*% beta`, which raises an error for a parameter of the wrong
length. It does not work for a vectorized expression, since R recycles a
shorter vector silently whenever its length divides, so
`sum((p - c(1, 2, 3))^2)` is a finite function of one parameter as well
as of three. The probe then signals an error that names the lengths
found, and `npar` has to be given.

## Bounds by reparametrization

Each bounded coordinate is reparametrized onto the whole real line, by a
shifted log for a one-sided bound and a scaled logit for a two-sided
one, and the optimizer runs unconstrained in the new variable. Every
point that it proposes is admissible by construction, so there is no
rejection step and no boundary for a line search to cross, and a method
needs no code of its own for the bounds.
[`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)
is the map.

`lower` and `upper` are two vectors, as in
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) and
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), so `lower = 0`
applies the lower bound zero to every parameter without one value per
coefficient. A length that is neither 1 nor `p` is rejected, and so is a
`lower` that is at or above its `upper`. The starting value must lie
strictly inside its bounds: `par = 0` with `lower = 0` is an error,
because the transformed coordinate is infinite there.

## An optimum lying on a bound

Reaching a bound exactly requires the transformed variable to run to
infinity, so a solution on a bound is approached and not attained. The
stopping rule ends the run, because the chain-rule factor \\\partial x /
\partial \eta\\ decays to zero as the bound is approached and takes the
transformed gradient with it. The run therefore reports
`converged = TRUE` at a point near the bound, and how near is set by the
tolerance rather than by the budget.

The bounds examples below show this on \\\sum (x - (1, 2))^2\\ with
`upper = c(5, 1)`, whose minimum over the box lies on the upper bound of
the second coordinate: a smaller tolerance of the rule gives a point
closer to the bound, after more iterations.

For the parameters for which the bounds are intended, such as a positive
variance or a probability inside the unit interval, the optimum is
usually interior and this does not arise. For a box-constrained problem
whose constraints are active at the solution, an active-set method is
more suitable.

## See also

[`maximize()`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
for the other direction,
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
for the algorithms,
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
for the stopping rules,
[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
for the starters,
[`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)
for the map that removes the bounds.

## Examples

``` r
q <- function(p) sum((p - c(1, 2))^2)
qg <- function(p) 2 * (p - c(1, 2))

# The first run supplies the gradient, and the second differences it.
# The two runs reach the same point and differ in the evaluation counts.
a <- minimize(bfgs(), q, par = c(0, 0), gr = qg)
b <- minimize(bfgs(), q, c(0, 0))
all.equal(a@par, b@par, tolerance = 1e-6)
#> [1] TRUE
rbind(supplied = unlist(a@counts), differenced = unlist(b@counts))
#>              f g h
#> supplied     3 3 0
#> differenced 15 0 0

# One bound applies to every parameter: a scale that must stay positive.
minimize(bfgs(), q, c(0.5, 0.5), lower = 0)@par
#> [1] 1 2

# Bounds may also be given per parameter. The unconstrained minimum is at (1, 2), so the
# second coordinate approaches its upper bound of 1 and stops just short of
# it, at a distance set by the tolerance.
minimize(bfgs(), q, c(0.5, 0.5), lower = c(0, 0), upper = c(5, 1))@par
#> [1] 1.0000005 0.9999996
minimize(bfgs(criterion = crit_grad(1e-12)), q, c(0.5, 0.5),
         lower = c(0, 0), upper = c(5, 1))@par
#> [1] 1 1

# Without a starting value, the bounds give the number of parameters.
minimize(bfgs(), q, start_zeros(), lower = c(0, 0), upper = c(5, 10))@par
#> [1] 1 2

# A start on its own bound is rejected: the transform is infinite there.
try(minimize(bfgs(), q, c(0, 0), lower = 0))
#> Error : The starting value for parameter 1 must lie strictly inside its bounds (0, Inf); it is 0.

# A gradient that belongs to a different model gives a warning, before the
# run, that names both rates, and the run then proceeds.
wrong <- minimize(gd(maxit = 2), q, c(0, 0), gr = function(p) 5 * qg(p))
#> Warning: 'gr' does not appear to be the gradient of 'fn': along the gradient direction at 'par',
#>   'fn' changes at rate 4.47 where 'gr' predicts 22.4. Check that the two
#>   compute the same model. options(optimizers7.check_gradient = FALSE) turns this check off.
wrong@converged
#> [1] FALSE
```
