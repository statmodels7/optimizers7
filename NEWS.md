# optimizers7 0.5.0

* `newton()` scales a direction that came out of a REPAIRED Hessian to
  `min(1, 1/||d||_inf)`, a displacement of order one in the parameters. Where
  the Cholesky succeeds nothing is touched: that step is the Newton step and
  its length is the objective's own curvature. Where it fails, the component
  in the floored subspace is `g_i / (floor * max|lambda|)`, whose size comes
  from the `floor` argument and from no curvature at all, so `step = 1` claims
  a Newton unit it does not have. The rule is the one `bb()` applies when its
  secant pair reports no curvature and the one `bfgs()` and `lbfgs()` apply to
  their first direction, and it ONLY EVER SHORTENS.

* The same scaling reaches the gradient `newton()` falls back on when the
  Hessian is not finite or a solve fails, which had never been given it.

* The trace reports a capped step as `"hessian modified (capped)"`, plain
  `"hessian modified"` where the direction was already of order one.

* Measured. On the package's own indefinite-Hessian test the repaired
  direction's first component is `5.03e+06`, and the line search was the only
  thing bounding it: the saddle costs 28 function evaluations before and 5
  after with `hessian_mod = "eigen"`, 30 and 5 with `"ridge"`. Over
  `test_problems()` under `newton()`, 743 evaluations become 480 and 8 of 8
  problems converge where 7 did -- `abs_sum` went from 180 evaluations and a
  gradient of 1.0 to 27 evaluations and a gradient of exactly 0. The four
  problems where the cap never fires (sphere, rosenbrock, booth, powell) are
  unchanged evaluation for evaluation, as are both Rosenbrock controls, whose
  Hessian is never repaired.

# optimizers7 0.4.0

* `armijo()`, `wolfe()` and `nonmonotone()` take a `resolution`: the smallest
  difference in the objective that means anything, in the objective's own
  units. It defaults to `0`, which does not ask the question, so nothing
  changes for a caller who does not set one.

* It may be a FUNCTION of no arguments rather than a number, for an objective
  whose resolution moves as the run goes -- a fit warm-started from the
  previous evaluation locates its own answer better each time, so the reading
  at the start is the one from the worst point of the whole run. Measured on a
  penalized smooth, the cold-start reading is 14 to 22 times the best later
  one. It is asked once per invocation of the search rather than once per
  trial, so it costs one call an iteration, and an answer that is not finite
  and positive is read as `0`.

* It is for an objective computed by a PROCEDURE rather than by a formula --
  a fit warm-started from wherever the last evaluation ended, a quadrature
  whose panels move, a simulation -- which returns slightly different values
  for the same argument. Below that spread its values carry no information,
  and a search that keeps backtracking spends its whole budget arriving back
  where it started: measured on a marginal likelihood whose coefficients are
  refitted at every evaluation, 30 evaluations, the backtracking budget
  exactly, each one a whole inner fit.

* **The question is asked once, of the full step, and not inside the
  backtracking loop.** The quantity is the improvement the method's own linear
  model predicts, `s0 * |g'd|`, and where that is below the resolution the
  search returns immediately: the point is optimal to the accuracy the
  objective has, reported as `no decrease above the objective's resolution`
  rather than as a stopping rule being met.

* ⚠️ Asking it inside the loop was tried first and is UNSAFE. There the two
  situations cannot be told apart, `x + s d` tending to `x` as the step shrinks
  whether the point is optimal or the DIRECTION is wrong; measured, a
  mis-stated gradient at a point nowhere near stationary was promoted to a
  converged run. Tested at the full step the two separate, a bad direction
  predicting a large improvement and still being reported as the failure it is.
  Both cases are in the tests.

* The quantity is the predicted decrease and not the Armijo demand
  `c1 * s0 * |g'd|`, which is four orders smaller and would fire where the
  method still had real progress to make.

# optimizers7 0.3.0

* `bfgs()` and `lbfgs()` scale their FIRST direction to a step of order one
  in the parameters. A quasi-Newton direction is scaled so that `step = 1`
  means the Newton step, which is why both default to it; on the first
  iteration -- and after a reset -- there is no curvature information and the
  direction degenerates to `-g`, for which one is not a natural unit at all.
  The trial displacement is then the gradient itself, so on a badly scaled
  objective the first point tried is an arbitrary distance away and the line
  search pays to backtrack all the way in.

  It is the rule `bb()` already uses when its secant pair reports no
  curvature, adopted there after being measured against three alternatives,
  and for the same two reasons: `1/||g||_inf` cannot freeze, not depending on
  the step it replaces, and cannot explode, scaling with the gradient. It
  only ever SHORTENS -- `min(1, 1/||g||_inf)` -- so a problem whose gradient
  at the start is already of order one takes the identical path, which is
  what makes this safe to change under the packages that depend on the
  default.

  `gd()` and `cg()` are deliberately unchanged: there the direction is the
  gradient by design and the step length is the line search's whole job, so
  `step = 1` was never claiming to be a Newton unit.

  Measured on `test_problems()`, before against after: the same 14 of 16
  converged, total evaluations 984 to 883, and the gradient at the reported
  solution equal or better on every problem (rosenbrock under bfgs, 1.15e-07
  to 4.37e-10). Per problem, rosenbrock 108 to 90, booth 26 to 16, powell 112
  to 84. On `0.5*1e4*|x|^2` from `x0 = 1`, where the arithmetic is explicit,
  19 evaluations become 4.

  What prompted it was statmodels7's marginal criterion over a score-driven
  panel, whose derivative grows with the number of penalized coordinates: a
  gradient of 12.3 on a log-scale hyperparameter put the first trial value at
  `exp(-12.3)` of the start, and the search spent itself backtracking through
  a region where the inner fit could not be evaluated. Four panels that
  needed a derivative-free search, or failed outright, now fit in 5 to 8
  criterion evaluations against 23 to 55.

# optimizers7 0.2.0

* `sa()`, simulated annealing with an adaptive step. The parameters move one
  coordinate at a time and each coordinate's step is adjusted every few
  sweeps to hold its acceptance rate near a target (Corana et al. 1987),
  which is what makes the method usable on a statistical objective, whose
  unconstrained coordinates sit on scales orders of magnitude apart. The
  proposal is uniform or Cauchy, the second being fast simulated annealing
  (Szu and Hartley 1987) and the `q = 2` member of the Tsallis family; the
  general Tsallis visiting distribution is deliberately absent, its generator
  being one that could not be validated against anything already here.
  The initial temperature is calibrated from the objective's own variation
  unless the caller sets it, so the same problem scaled by a million is
  solved as well. What the run returns is the best point SEEN and not the
  last, and `converged` is never inferred from the schedule finishing: the
  stationarity reported is Corana's own termination rule, so
  `crit_stationary()` is that rule rather than a second convention beside it.
  `sa_run_r()` is the R twin the compiled loop is held to -- the two draw
  from R's generator in the same order, so from one seed they are the same
  run and the test needs no tolerance. Measured, the port is worth 1.74x on
  an objective costing 0.7 microseconds and 1.34x on one costing 3.8; on an
  objective of half a millisecond, which is what a modelling layer's inner
  one costs, the loop's overhead is a quarter of a per cent and the port
  buys nothing.

* `chain()`, which runs optimizers one after another, each starting where the
  previous finished: `chain(sa(), lbfgs())` explores and then descends, and
  neither method knows about the other. It is the second wrapper of this
  shape after `multistart()`, and the two compose. Each stage carries its own
  criterion and budgets; the point, the value and `converged` are the last
  stage's, the work is summed, and the trace carries a `stage` column.

* `check_optimizer()`'s ninth check runs its `maximize` comparison from the
  SAME random stream as the `minimize` it is compared against, as the eighth
  check already re-seeded. Without it the check silently required the
  optimizer to converge tightly enough that two independent runs agree to
  1e-5, which the mesh-shrinking methods happen to do -- so it passed for the
  wrong reason everywhere and failed for the wrong reason on a global search.

# optimizers7 0.1.0

## Methods

* `prox_grad()`, the proximal gradient method: FISTA with backtracking and
  adaptive restart. It is the one algorithm written in R rather than C++,
  every iteration calling a caller-supplied proximal operator. Acceleration
  matters on ill-conditioned problems and not otherwise: at a condition
  number of 3 the plain iteration wins narrowly, at 55 it takes 4153
  iterations against 126, and at 480 it does not converge in 50000 while
  the accelerated one takes 334.

* `optimizer_bounded()` declares whether a method accepts box bounds beside
  the objective. It is `FALSE` for `prox_grad()`, which takes its
  constraint inside the operator, where it composes with the term already
  there; `check_optimizer()` consults it before testing bounds.

* `nonmonotone()`, a third `line_search` subclass implementing the
  condition of Grippo, Lampariello and Lucidi, and the default for `bb()`.
  With `memory = 0` it is `armijo()` value for value: on Rosenbrock, 68
  iterations and 77 evaluations against 82 and 186, with eleven uphill
  steps accepted out of sixty-seven where Armijo accepts none.

* `bb()` takes its fallback step from `1/max|g|` when the secant pair
  reports no usable curvature. Both constants one reaches for first are
  absorbing states: freezing the previous alpha traps a short step in its
  own shortness (873 of 945 iterations, on bb2), and restarting at `alpha0`
  does the same wherever `alpha0` is itself too short (1395 rejected steps
  in 1521 iterations, on a boxed quadratic). `alpha_max`, which the SPG
  literature uses, hands the search a direction of length 1e10 and the run
  stopped a unit from the solution reporting success. The curvature test is
  also relative now, `s'y > curv_tol*||s||*||y||`, as `bfgs()` already had
  it.

* `start_zeros()` and `start_runif()`, so that a starting value need not be
  written out, and `multistart()` runs its starts in parallel by default.

## Stopping and reporting

* The default gradient tolerance is `1e-6`, in `crit_grad()` and in every
  method that restates a criterion. A line search accepts a step only when
  the objective decreases by a definite amount, so the smallest gradient a
  run can reach is about `sqrt(2*lambda*eps*|f*|)`, which grows with the
  value at the solution. Every problem in `test_problems()` has `f* = 0`,
  where that floor is itself zero; a log-likelihood is of order one at its
  maximum and the floor is then around 1e-8. Measured across families,
  methods and seeds it reaches 1.06e-8.

* `gd()`, `cg()`, `bb()`, `newton()`, `bfgs()` and `lbfgs()` no longer
  restate the tolerances in their defaults, so a change to `crit_grad()`
  reaches them.

* A run asks its stopping rule before reporting a line-search failure. A
  search that could find no acceptable step broke the loop with
  `converged = FALSE` even when the start was already the answer, which is
  now the ordinary case for a caller handing in a closed-form estimate. The
  criterion is asked with `have_old = false`, so a rule reading a change in
  the objective returns `FALSE` by construction and only the state at the
  point can end the run.

* The gradient check is skipped at a stationary start. Differencing the
  objective along a direction of no slope gives its own truncation error,
  so the check warned at precisely the caller who supplied the exact
  optimum.

* A supplied gradient is checked once against the objective, one central
  difference along the gradient direction at `par`, and a gross
  disagreement draws a warning naming both rates. Disable with
  `options(optimizers7.check_gradient = FALSE)`.

* `newton()` warns when the budget admits fewer than two iterations, and
  `max_eval` defaults to `Inf`.

## Interface

* Box constraints are two vectors, `lower` and `upper`, rather than a list
  of pairs, so `lower = 0` says that every parameter is positive.

* `adam(resample = )` and `finite_sum()` are removed. An objective that
  draws its own minibatches is a closure the caller writes in one line, and
  between them they cost a second objective class, a rule for which
  criteria they permitted, a token in the criterion machinery, a branch in
  the compiled loop and an argument to select the path.

* `cpp_objective()` is removed. The pointer type
  `double(*)(const arma::vec&)` has no closure, so any real objective
  needed C++ globals, which is not reentrant.

* `multistart()` runs its starts sequentially under `pkgload`. Its PSOCK
  workers load the installed copy of the package, so an S7 optimizer built
  in a development namespace dispatched against the installed namespace's
  methods and the inner runs came back wrong without an error.

## Validation and documentation

* `check_optimizer()` separates what an optimizer promises from how well it
  does it: twelve checks on what it reports -- that `value` is the
  objective at `par`, that `converged` follows the rule and is never
  inferred from the run ending, that bounds hold strictly -- and then a
  battery of standard problems whose gaps are printed as information. It
  found two defects in `bundle()` the afternoon it was written, one of them
  a stationarity measure containing the trust parameter, which the
  algorithm is free to shrink for reasons of its own.

* `check_bounds()`, `check_criterion()` and `as_objective()` are exported,
  because writing the extension vignette showed that a user could not
  otherwise write a conforming optimizer: honoring the box constraints and
  rejecting a stopping rule the method cannot evaluate both needed internal
  functions.

* Every constructor cites the paper its method comes from, every exported
  topic carries a `\value` and an executable example, and the class pages
  live at `X-class` with the plain name as an alias, two topics differing
  only in case being one file on Windows.

* A vignette on extending the package, a README with badges, a pkgdown site
  and continuous integration on five platforms.
