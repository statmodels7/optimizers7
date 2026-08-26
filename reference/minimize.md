# Minimize a Function

Runs an optimizer on an objective and returns the point it reached, the
value there, and the rule that stopped it. This is the one entry point:
every algorithm in the package is a method of this generic, so changing
method means changing the first argument and nothing else. Box
constraints are accepted by every method. Everything here minimizes;
[`maximize()`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
negates the objective and restores the sign of the answer.

## Usage

``` r
minimize(
  optimizer,
  fn,
  par,
  gr = NULL,
  he = NULL,
  lower = -Inf,
  upper = Inf,
  ...
)
```

## Arguments

- optimizer:

  An
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  object carrying the algorithm and its settings. The only argument
  dispatch reads.

- fn:

  The objective: a function of the parameter vector returning a single
  number, to be minimized. An object of another class works whenever
  [`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md)
  has a method for it.

- par:

  A numeric vector of starting values, or a starter object; see Starters
  below. With bounds it must lie **strictly** inside them.

- gr:

  The gradient, a function of the parameter vector. `NULL`, the default,
  has it differenced from `fn`. Ignored when `fn` is an objective
  carrying its own gradient.

- he:

  The Hessian, a function of the parameter vector. `NULL` is the default
  and only
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
  reads one; every other method accepts it and ignores it, so calling
  code need not branch on the algorithm.

- lower, upper:

  Box constraints, numeric of length 1 (applying to every parameter) or
  of length `p` (one each). Any other length is an error naming both.
  The defaults `-Inf` and `Inf` are no constraint, and a coordinate
  whose pair is infinite on both sides is left untransformed. `lower`
  must be strictly below `upper` in every coordinate.

- ...:

  Passed to the method dispatched on. No shipped method reads anything
  from it.

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
it stopped together with the rule that stopped it. Convergence is what a
stopping rule confirmed and is never inferred from the run having ended,
so a run that exhausts its budget comes back with `converged = FALSE`.

Dispatch is on `optimizer` alone, so each algorithm is written once. The
objective is normalized separately, by
[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md),
which dispatches on `fn`; a caller with its own kind of objective
registers one method there and every algorithm accepts it.

## Derivatives, supplied and differenced

A gradient that is not supplied is computed by central finite
differences, at a cost of \\2p\\ evaluations of the objective each time.
The result's `message` records that, so a run is never silently less
exact than it looks.

A gradient that *is* supplied is checked once, before the run: one
central difference along the gradient direction at `par`, two
evaluations. A gross disagreement draws a warning naming both rates, as
in

    'gr' does not appear to be the gradient of 'fn': along the gradient
    direction at 'par', 'fn' changes at rate 4.47 where 'gr' predicts 22.4.

and the run proceeds. A `gr` computed from a different model than `fn`
is otherwise very hard to see: it surfaces as a mute line-search failure
at the first iteration. Set
`options(optimizers7.check_gradient = FALSE)` to turn the check off.

## Starters

`par` may be a **starter** object instead of a vector:
[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
for all zeros,
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md)
for a uniform draw from a chosen range. Both work on the unconstrained
scale and are mapped back through the bounds, so one constant means
something sensible for every kind of parameter (zero becomes one for a
variance, one half for a probability) and no draw can land outside its
box.

A starter has to be told how many parameters there are, in one of three
ways. Say so with `start_zeros(npar = 3)`. Failing that, a `lower` or
`upper` of more than one element answers the question, bounds being one
per parameter. Failing that, the objective is probed by
[`infer_npar()`](https://statmodels7.github.io/optimizers7/reference/infer_npar.md),
which tries lengths 1 to 50 until one is accepted, once, before the run
begins.

That last route settles any objective with a fixed width built into it,
which is most real ones: `X %*% beta` with a parameter of the wrong
length is an error. It cannot settle a vectorized toy, since R recycles
a shorter vector silently whenever its length divides, so
`sum((p - c(1, 2, 3))^2)` is a perfectly finite function of one
parameter as well as of three. It then rejects, naming the lengths it
found, and asks for `npar`.

## Bounds are removed, not enforced

Each bounded coordinate is reparametrized onto the whole real line, by a
shifted log for a one-sided bound and a scaled logit for a two-sided
one, and the optimizer runs unconstrained in the new variable. Every
point it proposes is admissible by construction, so there is no
rejection step and no boundary for a line search to trip over, and any
method takes bounds without knowing they exist.
[`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)
is the map.

`lower` and `upper` are two vectors, as in
[`stats::optim()`](https://rdrr.io/r/stats/optim.html) and
[`stats::nlminb()`](https://rdrr.io/r/stats/nlminb.html), so `lower = 0`
says that every parameter is positive without writing out one pair per
coefficient. A length that is neither 1 nor `p` is refused, and so is a
`lower` at or above its `upper`. The **starting value must lie strictly
inside its bounds**: `par = 0` with `lower = 0` is an error, the
transformed coordinate being infinite there.

## An optimum lying on a bound

Reaching a bound exactly requires the transformed variable to run to
infinity, so a solution on a bound is approached and not attained. What
stops the run is the stopping rule, because the chain-rule factor
\\\partial x / \partial \eta\\ decays to zero as the bound is approached
and takes the transformed gradient with it. The run therefore reports
`converged = TRUE` at a point near the bound, and how near is set by the
tolerance rather than by the budget.

Measured on \\\sum (x - (1, 2))^2\\ with `upper = c(5, 1)`, whose second
coordinate wants to sit at its ceiling: the gap to the bound is
`3.6e-07` at `crit_grad(1e-6)`, `1.1e-11` at `1e-10` and `2.7e-15` at
`1e-14`, and raising `maxit` from 50 to 5000 changes nothing, all three
runs stopping after 22 iterations.

For the statistical use this exists to serve, a positive variance or a
probability inside the unit interval, the optimum is interior and none
of this arises. For a genuine box-constrained problem whose constraints
are active at the solution, an active-set method is the right tool.

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

# With the gradient supplied, and without, so it is differenced. The two
# reach the same point; only the evaluation counts differ.
a <- minimize(bfgs(), q, par = c(0, 0), gr = qg)
b <- minimize(bfgs(), q, c(0, 0))
all.equal(a@par, b@par, tolerance = 1e-6)
#> [1] TRUE
rbind(supplied = unlist(a@counts), differenced = unlist(b@counts))
#>              f g h
#> supplied     3 3 0
#> differenced 15 0 0

# One bound for every parameter: a scale that must stay positive.
minimize(bfgs(), q, c(0.5, 0.5), lower = 0)@par
#> [1] 1 2

# Or one per parameter. The unconstrained minimum is at (1, 2), so the
# second coordinate is pushed against its ceiling of 1 and stops just
# short of it, by an amount the tolerance sets.
minimize(bfgs(), q, c(0.5, 0.5), lower = c(0, 0), upper = c(5, 1))@par
#> [1] 1.0000005 0.9999996
minimize(bfgs(criterion = crit_grad(1e-12)), q, c(0.5, 0.5),
         lower = c(0, 0), upper = c(5, 1))@par
#> [1] 1 1

# No starting value at all: the bounds say there are two parameters.
minimize(bfgs(), q, start_zeros(), lower = c(0, 0), upper = c(5, 10))@par
#> [1] 1 2

# A start sitting on its own bound is refused: the transform is infinite
# there.
try(minimize(bfgs(), q, c(0, 0), lower = 0))
#> Error : The starting value for parameter 1 must lie strictly inside its bounds (0, Inf); it is 0.

# A gradient belonging to a different model draws a warning naming both
# rates, before the run, and the run then proceeds.
wrong <- minimize(gd(maxit = 2), q, c(0, 0), gr = function(p) 5 * qg(p))
#> Warning: 'gr' does not appear to be the gradient of 'fn': along the gradient direction at 'par',
#>   'fn' changes at rate 4.47 where 'gr' predicts 22.4. Check that the two
#>   compute the same model. options(optimizers7.check_gradient = FALSE) turns this check off.
wrong@converged
#> [1] FALSE
```
