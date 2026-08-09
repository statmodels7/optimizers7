# Changelog

## optimizers7 0.0.0.9000

### Methods

- [`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md),
  the proximal gradient method: FISTA with backtracking and adaptive
  restart. It is the one algorithm written in R rather than C++, every
  iteration calling a caller-supplied proximal operator. Acceleration
  matters on ill-conditioned problems and not otherwise: at a condition
  number of 3 the plain iteration wins narrowly, at 55 it takes 4153
  iterations against 126, and at 480 it does not converge in 50000 while
  the accelerated one takes 334.

- [`optimizer_bounded()`](https://statmodels7.github.io/optimizers7/reference/optimizer_bounded.md)
  declares whether a method accepts box bounds beside the objective. It
  is `FALSE` for
  [`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md),
  which takes its constraint inside the operator, where it composes with
  the term already there;
  [`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
  consults it before testing bounds.

- [`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md),
  a third `line_search` subclass implementing the condition of Grippo,
  Lampariello and Lucidi, and the default for
  [`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md).
  With `memory = 0` it is
  [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
  value for value: on Rosenbrock, 68 iterations and 77 evaluations
  against 82 and 186, with eleven uphill steps accepted out of
  sixty-seven where Armijo accepts none.

- [`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md)
  takes its fallback step from `1/max|g|` when the secant pair reports
  no usable curvature. Both constants one reaches for first are
  absorbing states: freezing the previous alpha traps a short step in
  its own shortness (873 of 945 iterations, on bb2), and restarting at
  `alpha0` does the same wherever `alpha0` is itself too short (1395
  rejected steps in 1521 iterations, on a boxed quadratic). `alpha_max`,
  which the SPG literature uses, hands the search a direction of length
  1e10 and the run stopped a unit from the solution reporting success.
  The curvature test is also relative now, `s'y > curv_tol*||s||*||y||`,
  as
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  already had it.

- [`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
  and
  [`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md),
  so that a starting value need not be written out, and
  [`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
  runs its starts in parallel by default.

### Stopping and reporting

- The default gradient tolerance is `1e-6`, in
  [`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
  and in every method that restates a criterion. A line search accepts a
  step only when the objective decreases by a definite amount, so the
  smallest gradient a run can reach is about `sqrt(2*lambda*eps*|f*|)`,
  which grows with the value at the solution. Every problem in
  [`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
  has `f* = 0`, where that floor is itself zero; a log-likelihood is of
  order one at its maximum and the floor is then around 1e-8. Measured
  across families, methods and seeds it reaches 1.06e-8.

- [`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md),
  [`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md),
  [`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md),
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  and
  [`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
  no longer restate the tolerances in their defaults, so a change to
  [`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
  reaches them.

- A run asks its stopping rule before reporting a line-search failure. A
  search that could find no acceptable step broke the loop with
  `converged = FALSE` even when the start was already the answer, which
  is now the ordinary case for a caller handing in a closed-form
  estimate. The criterion is asked with `have_old = false`, so a rule
  reading a change in the objective returns `FALSE` by construction and
  only the state at the point can end the run.

- The gradient check is skipped at a stationary start. Differencing the
  objective along a direction of no slope gives its own truncation
  error, so the check warned at precisely the caller who supplied the
  exact optimum.

- A supplied gradient is checked once against the objective, one central
  difference along the gradient direction at `par`, and a gross
  disagreement draws a warning naming both rates. Disable with
  `options(optimizers7.check_gradient = FALSE)`.

- [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
  warns when the budget admits fewer than two iterations, and `max_eval`
  defaults to `Inf`.

### Interface

- Box constraints are two vectors, `lower` and `upper`, rather than a
  list of pairs, so `lower = 0` says that every parameter is positive.

- `adam(resample = )` and `finite_sum()` are removed. An objective that
  draws its own minibatches is a closure the caller writes in one line,
  and between them they cost a second objective class, a rule for which
  criteria they permitted, a token in the criterion machinery, a branch
  in the compiled loop and an argument to select the path.

- `cpp_objective()` is removed. The pointer type
  `double(*)(const arma::vec&)` has no closure, so any real objective
  needed C++ globals, which is not reentrant.

- [`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
  runs its starts sequentially under `pkgload`. Its PSOCK workers load
  the installed copy of the package, so an S7 optimizer built in a
  development namespace dispatched against the installed namespace’s
  methods and the inner runs came back wrong without an error.

### Validation and documentation

- [`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
  separates what an optimizer promises from how well it does it: twelve
  checks on what it reports – that `value` is the objective at `par`,
  that `converged` follows the rule and is never inferred from the run
  ending, that bounds hold strictly – and then a battery of standard
  problems whose gaps are printed as information. It found two defects
  in
  [`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
  the afternoon it was written, one of them a stationarity measure
  containing the trust parameter, which the algorithm is free to shrink
  for reasons of its own.

- [`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md),
  [`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
  and
  [`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md)
  are exported, because writing the extension vignette showed that a
  user could not otherwise write a conforming optimizer: honoring the
  box constraints and rejecting a stopping rule the method cannot
  evaluate both needed internal functions.

- Every constructor cites the paper its method comes from, every
  exported topic carries a `\value` and an executable example, and the
  class pages live at `X-class` with the plain name as an alias, two
  topics differing only in case being one file on Windows.

- A vignette on extending the package, a README with badges, a pkgdown
  site and continuous integration on five platforms.
