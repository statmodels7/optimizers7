# optimizers7 0.15.0

## Repairs

* In the methods with a line search (`newton()`, `chebyshev()`, `bfgs()`,
  `lbfgs()`, `gd()`, `cg()`, `bb()`), an accepted step that moves the
  iterate by no more than `4 * eps * max(1, max(abs(x)))` is treated as a
  line search that found no acceptable step: the stopping rule is evaluated
  without the previous iterate, and otherwise the run stops with
  `converged = FALSE`, the message `the step did not move the iterate
  beyond rounding` and the safeguard `step below rounding`. In the same
  methods, when the stopping rule is met on a step `dx` with
  `abs(sum(g * dx)) > max(1, abs(f))` at the new point, the step is undone
  and handled in the same way, with the message `the step that met the
  stopping rule is not resolved by the objective` and the safeguard `step
  not resolved`. Where the objective is rounding noise, such steps passed
  the Armijo condition when the noise went down, and `crit_abs_par()` or
  `crit_abs_obj()` then reported convergence at a gradient of order 1e27
  (the exploded score filter of a beta share fitted by statmodels7, where
  the flag differed between platforms at the same point: 30 to 34 of 200
  perturbed starts converged, none now).

# optimizers7 0.14.0

## Repairs

* `brent()` reports its own stop under its own label. When the bracket has
  shrunk to the resolution of the iterate, the run ends with
  `converged = TRUE` and `criterion_met` equal to `bracket at machine
  resolution`; it previously reported the label of the stopping rule, even
  for `crit_never()`. The state passed to the rule now has `iter` and
  `x_old`, so `crit_abs_par()` and `crit_rel_par()` can fire. When
  `max_expand` expansions or the evaluation budget end the bracketing phase
  before a minimum is bracketed, the run stops at the lowest point evaluated
  with `converged = FALSE` and a message; it previously reduced a bracket
  that contained no minimum and reported convergence (9.47 on
  `(x - 1000)^2` with `max_expand = 3`).
* `prox_grad()` counts the evaluations of a differenced gradient in
  `counts` and in `max_eval` (54 reported against 446 calls before). The
  trace has the columns `iteration`, `value`, `gnorm`, `step` and
  `safeguard`, the last recording `step shortened` and `momentum restart`,
  so `summary()` counts them. The state passed to the rule has `iter`. With
  `accelerate = FALSE` the reported gradient is the proximal gradient
  mapping at the reported point; it was the mapping at the previous iterate.
  The gradient read at the iterate is reused when the next step starts
  there, so the plain method evaluates one gradient per iteration.
* `adam()` tests its stopping rule at the point produced by the last
  iteration of the budget. A run with `maxit` equal to the iteration at
  which the rule holds now reports `converged = TRUE`, as a run with a
  larger budget does.
* `nelder_mead()` rebuilds a degenerate simplex with edges equal to the size
  that it had reached. The edges were multiplied again by the size of the
  coordinates, so a simplex far from the origin was rebuilt about
  \|x\| times too large.
* `bundle()` stores the aggregate linearization at the stability center at
  which it was formed. On a serious step followed by a compression it was
  stored at the new center, where it can lie above the objective: with
  `bundle_size = 2` the sphere stopped at 0.044 reporting convergence. On
  the test problems, with `bundle_size` between 2 and 20, the repaired
  method converges on three runs that failed before and fails on two that
  converged.
* `maximize()` negates the third-derivative contraction `t3` of
  `chebyshev()`, as it negates `gr` and `he`.
* `multistart()` answers `optimizer_bounded()` with the answer of its inner
  optimizer, so a multistart over `prox_grad()` no longer accepts bounds
  that the inner runs reject. `with_trace()` has a method for `chain()` that
  switches on the trace of every stage.
* The counts are checked to be finite: `lbfgs(memory = Inf)`, which ended
  the R session, `bfgs(max_skip = Inf)`, `armijo(max_step = Inf)`,
  `nonmonotone(memory = Inf)`, `sa(cycles = Inf)`, `bundle(bundle_size =
  Inf)` and `multistart(ncores = Inf)` are rejected.
* The error for an invalid tolerance names the argument that was supplied
  (`'spread'`, `'alpha0'`, ...) and not `'tol'`.
* `format_elapsed()` chooses the unit after rounding, so that it no longer
  prints `59 min 60 s` or `1 h 60 min`.
* The analytic Hessian of the `beale` test problem is finite at `y = 0`.
* The message of a multistart whose first start failed no longer ends the
  quoted error with two periods, and the error raised when every start
  fails reads `The first error was:`.
* The warning of `newton()` for a small evaluation budget states the cost
  as `4*p^2 + 2*p`, the quantity that it computes.

## Documentation

* The documentation, the README, the vignette and the reference index were
  rewritten in a plain register, and claims that the code contradicted were
  corrected. Among them: `chebyshev()` reads the Hessian and shares the
  default rule of the gradient methods; `brent()` provides `"stationarity"`;
  the size of the Nelder-Mead simplex is the largest distance from the best
  vertex; the Corana measure of `sa()` and the cost of its temperature
  calibration (`max(10, 10 * p)` evaluations); the `floor` of `newton()` and
  the initial poll size of `compass()`; the trace labels of `cg()` and
  `bb()`; the conditions under which `converged` is `TRUE`; the seeds that a
  result records; the condition numbers of `rosenbrock` and `beale`
  (2508 and 162).

# optimizers7 0.13.0

* Every problem of `test_problems()` carries its analytic Hessian as `he`,
  and `check_optimizer()` passes it to the method under test. A method that
  requires a Hessian could not be checked before, since the battery supplied
  only the objective and the gradient.

# optimizers7 0.12.1

* The check that a supplied gradient belongs to the objective repeats its
  directional difference at steps ten and a hundred times smaller before it
  warns, and warns only if the disagreement persists. The step is scaled by
  the length of the starting vector, so a large coordinate made it too long
  for a strongly curved direction: every `gas()` fit in statmodels7 gave the
  warning (a difference of -0.00106 against an exact 0.00074 on the Nile
  flows), with a gradient that agrees with `numDeriv` to seven digits.

# optimizers7 0.12.0

* `brent()`, Brent's method for a function of one variable: the minimum is
  bracketed by golden expansion from the starting value and the bracket is
  reduced by parabolic interpolation, with a golden-section step wherever
  the parabola is not trusted. The stationarity measure is the half-width of
  the bracket, so `crit_stationary()` reads a tolerance on the scale of the
  parameter. Bounds are removed by reparametrization, as in the other
  methods. On (x - 2)^2 + exp(-x) from 0 it takes 12 evaluations against the
  42 of `nelder_mead()`, and it rejects a starting value of more than one
  coordinate.

# optimizers7 0.11.0

* `chebyshev()`, Newton's method with Chebyshev's third-order correction.
  Where the Hessian is positive definite the Newton step is corrected by the
  third derivative of the objective contracted twice with it,
  `d = dN - H^-1 T[dN, dN] / 2`, read from `t3(x, d)` passed to `minimize()`,
  or from one second difference of the gradient along `dN` without it. A
  correction longer than `ratio` times the Newton step, or one that is not a
  descent direction, is rejected and the Newton step kept; the rest is as in
  `newton()`. On Rosenbrock the method takes 18, 4 and 5 iterations against
  the 21, 14 and 8 of `newton()` from three starts. A `t3` together with
  `lower` or `upper` is rejected, since the run is then on the unconstrained
  scale.

# optimizers7 0.10.0

* `newton()` takes `typical`, the typical size of each parameter (the scaled
  norm of Dennis and Schnabel). Every length used by the method is then
  measured in these units, the largest |d_i| / t_i in place of the largest
  |d_i|: the bound `max_length`, the length-one scaling of a direction whose
  Hessian was repaired, and the scaling of a gradient fallback. `NULL`, the
  default, measures them in the units of the parameters, and a run that does
  not set it is unchanged, which a test checks on Rosenbrock to the
  evaluation counts. On a restricted likelihood in a variance near 236,
  started at 651 where it is not concave, the repaired steps had moved the
  variance by one unit per iteration.

# optimizers7 0.9.0

* `newton()` takes `max_length`, the largest component that a step may
  have. A direction longer than this in the infinity norm is scaled down to
  this length before the line search, and the trace records the step as
  `newton step capped`. The default is `Inf`, which leaves every step as it
  is, so a run that does not set the argument is the run of the previous
  release.

  The bound is the `maxstep` of the line-search Newton method of Dennis and
  Schnabel (1983) and the `maxNstep` of the smoothing-parameter iteration of
  mgcv. Where the objective flattens towards an asymptote the curvature is
  small and the Newton step long, so the step can carry the iterate past the
  minimum into the flat region, where the line search accepts it and the
  gradient is close to zero. On the REML criterion of statmodels7 for a
  ridge penalty on the fifteen standardized predictors of `MASS::UScrime`,
  the first step from lambda = 1 had a length of 23.5 on the log scale and
  reached lambda = 1.6e10, where the criterion is flat at -65.40, while its
  maximum, -53.30, lies near lambda = 91.5. With `max_length = 5` the same
  search reaches the maximum.

  Near a minimum the Newton steps are short and the bound does not bind, so
  the quadratic convergence of the method is kept. On Rosenbrock from the
  customary start a bound of 10 leaves the run identical to the unbounded
  one, to the evaluation counts. On sqrt(1 + x^2), whose Newton step from 2
  is -10, a bound of 1 takes every step at full length and needs fewer
  evaluations than the unbounded method, which backtracks its first step.

# optimizers7 0.8.0

* The adaptive restart of `prox_grad()` compares an increase in the
  objective with an allowance for its rounding error, eight units in the
  last place of the current objective, and not with zero. The objective is
  a sum, so its error grows with the number of terms, and a bare `>` test
  fires on that error near the solution, discarding the momentum at each
  spurious reset.

  On the lasso of the `prox_grad()` page (200 observations, eight
  coefficients, three non-zero in the truth) at `crit_grad(1e-9)` the run
  took 21646 iterations before and 782 after, at the same support and the
  same objective. At the default `crit_grad(1e-6)` both take 11. Allowances
  between one and thirty-two units gave the same run, and 256 units cost
  iterations on an ill-conditioned problem by suppressing real restarts. On
  a smooth quadratic in eight unknowns at `crit_grad(1e-8)` the allowance
  changes the iteration counts of the restarted method from 149, 571 and
  934 to 150, 573 and 1042 at condition numbers 55, 480 and 2400.

  `restart = FALSE` is unchanged: 89 iterations on the same lasso.

# optimizers7 0.7.0

* The abstract `starter` class is exported. A stopping rule inherits from
  `criterion()` and an algorithm from `optimizer()`, both exported, while a
  starting-value generator had to inherit from an internal class:
  `starting_values()` was exported, but `minimize()` rejected the resulting
  object with `'par' must be a numeric vector of starting values`, since it
  tests inheritance from the class. A user-defined starter now runs as the
  two shipped ones do. Existing calls are unchanged.

* `bfgs()` and `lbfgs()` validate `curv_tol`, and `bfgs()` validates
  `max_skip`; both were accepted at any value. `curv_tol` multiplies
  `|s||y|` in the test of the curvature of a secant pair, so a negative value
  weakens the protection (`bfgs(curv_tol = -1)` ran identically to the
  default on Rosenbrock, the protection never firing). `max_skip` is
  compared with `>=`, so a negative value behaves as zero (`max_skip = -3`,
  `0` and `5` gave the same 500 iterations when the skip path was forced).
  Zero is admitted for both: `curv_tol = 0` is the textbook condition
  `s'y > 0`, and `max_skip = 0` resets the approximation on the first
  skipped update.

* `check_count()` rejects a fractional value, as its page stated. As a
  consequence `lbfgs(memory = 2.5)` is rejected instead of running at a
  memory of 2, and `armijo(max_step = 2.5)`, `wolfe(max_step = 2.5)` and
  `nonmonotone(max_step = 2.5)` are rejected, as the example of `armijo()`
  showed.

* `check_whole()` and `check_nonneg()` validate a count that may be zero and
  a threshold that may not be negative. `nonmonotone()` calls the first in
  place of its inline check, with the same message and bounds.

# optimizers7 0.6.0

* The six gradient methods `gd()`, `cg()`, `bb()`, `newton()`, `bfgs()` and
  `lbfgs()` default to `crit_any(crit_grad(), crit_abs_obj(),
  crit_abs_par())` instead of `crit_any(crit_grad(), crit_rel_obj())`. A
  disjunction only becomes weaker as terms are added, so a run that stops
  under the new rule would not have stopped earlier under the old one. Over
  `test_problems()`, six methods on eight problems, 44 of 48 runs converge
  instead of 41, at 19370 evaluations of the objective and the gradient
  instead of 22299. The three runs gained are `cg` and `bb` on the
  non-smooth `abs_sum` and `gd` on Beale; none is lost.

  A run stops sooner and reports a point further from the solution. Ten of
  the 48 end at a gradient more than a hundred times larger: `bfgs` on
  Rosenbrock ends at 8.1e-06 instead of 4.4e-10, the objective 3.2e-13 above
  its minimum instead of 1.2e-21. One run differs in substance: `cg` on
  `abs_sum` ends 4.5e-02 above the minimum instead of 1.9e-03 and now
  reports convergence, since on a non-smooth problem a smooth method stalls
  far from the solution and a rule that reads a stall reports it as
  convergence. `criterion = crit_grad()` restores the previous behavior.

* `crit_rel_obj()` leaves the default rule, since it never fired in it: over
  the same 48 runs the rule with it and the rule without it agree on every
  convergence flag, every evaluation count and every reported point. It
  remains exported.

# optimizers7 0.5.0

* `newton()` multiplies a direction computed from a repaired Hessian by
  `min(1, 1/||d||_inf)`, which gives a displacement of order one in the
  parameters. Where the Cholesky factorization succeeds the direction is
  unchanged, since it is the Newton step and its length is set by the
  curvature. Where it fails, the component in the floored subspace is
  `g_i / (floor * max|lambda|)`, whose size comes from `floor` and not from
  the curvature, so `step = 1` is not a Newton unit there. The rule is the
  one that `bb()` applies when its secant pair has no usable curvature and
  that `bfgs()` and `lbfgs()` apply to their first direction, and it never
  lengthens a step.

* The same scaling applies to the gradient direction that `newton()` uses
  when the Hessian is not finite or a solve fails.

* The trace reports a scaled step as `"hessian modified (capped)"`, and
  `"hessian modified"` where the direction was already of order one.

* On the indefinite-Hessian test of the package the first component of the
  repaired direction was `5.03e+06`, bounded only by the line search: the
  saddle cost 28 evaluations before and 5 after with `hessian_mod = "eigen"`,
  30 and 5 with `"ridge"`. Over `test_problems()` under `newton()`, 743
  evaluations became 480 and 8 of 8 problems converged instead of 7
  (`abs_sum` went from 180 evaluations and a gradient of 1.0 to 27
  evaluations and a gradient of exactly 0). The four problems on which the
  scaling never applies (sphere, rosenbrock, booth, powell) are unchanged to
  the evaluation counts.

# optimizers7 0.4.0

* `armijo()`, `wolfe()` and `nonmonotone()` take a `resolution`: the
  smallest difference in the objective that the objective can resolve, in
  its own units. It defaults to `0`, which disables the test, so a caller
  who does not set it sees no change.

* It may be a function of no arguments instead of a number, for an
  objective whose resolution changes during the run (a fit warm-started from
  the previous evaluation becomes more precise at each one, so the value
  measured at the start is the worst of the run; on a penalized smooth it was
  14 to 22 times the best later one). The function is called once per
  invocation of the search, so it costs one call per iteration, and a value
  that is not finite and positive is read as `0`.

* It is for an objective computed by a procedure (a warm-started fit, a
  quadrature whose panels move, a simulation), which returns slightly
  different values for the same argument. Differences below that spread
  carry no information, and a search that keeps backtracking uses its whole
  budget: on a marginal likelihood whose coefficients are refitted at every
  evaluation, 30 evaluations, the backtracking budget, each a complete inner
  fit.

* The test is made once, on the full step, and not inside the backtracking
  loop. The quantity is the decrease predicted by the linear model of the
  method, `alpha0 * |g'd|`, and when it is below the resolution the search
  returns at once: the point is optimal to the accuracy of the objective,
  and the result reports `no decrease above the objective's resolution`
  instead of a stopping rule being met.

* Inside the loop the test is unsafe, because the two situations cannot be
  distinguished as the step shrinks: `x + alpha d` tends to `x` whether the
  point is optimal or the direction is wrong, and an incorrect gradient at a
  point far from stationarity was reported as convergence. At the full step
  they differ, since a wrong direction predicts a large decrease and the
  search reports a failure. The tests cover both cases.

* The quantity is the predicted decrease and not the Armijo threshold
  `c1 * alpha0 * |g'd|`, which at the default `c1 = 1e-4` is four orders of
  magnitude smaller and would fire while the method can still make
  progress.

# optimizers7 0.3.0

* `bfgs()` and `lbfgs()` scale their first direction to a step of order one
  in the parameters. A quasi-Newton direction is scaled so that `step = 1`
  is the Newton step, which is why both default to it; on the first
  iteration, and after a reset, there is no curvature information and the
  direction is `-g`, for which one is not a natural unit. The trial
  displacement is then the gradient itself, so on a badly scaled objective
  the first point tried is at an arbitrary distance and the line search
  backtracks.

  The rule is the one that `bb()` uses when its secant pair has no
  curvature: `1/||g||_inf` cannot freeze, since it does not depend on the
  step that it replaces, and cannot explode, since it scales with the
  gradient. The factor `min(1, 1/||g||_inf)` never lengthens a step, so a
  problem whose gradient at the start is of order one takes the same path.

  `gd()` and `cg()` are unchanged: there the direction is the gradient by
  design and the step length is the task of the line search.

  On `test_problems()`, before against after: the same 14 of 16 converged,
  total evaluations 984 to 883, and the gradient at the reported solution
  equal or smaller on every problem (rosenbrock under bfgs, 1.15e-07 to
  4.37e-10). Per problem, rosenbrock 108 to 90, booth 26 to 16, powell 112
  to 84. On `0.5*1e4*|x|^2` from `x0 = 1`, 19 evaluations became 4. On the
  marginal criterion of statmodels7 over a score-driven panel, four panels
  that needed a derivative-free search, or failed, fit in 5 to 8 criterion
  evaluations against 23 to 55.

# optimizers7 0.2.0

* `sa()`, simulated annealing with an adaptive step. The parameters move one
  coordinate at a time and the step of each coordinate is adjusted every
  few sweeps to hold its acceptance rate near a target (Corana et al. 1987),
  which makes the method usable on a statistical objective, whose
  unconstrained coordinates lie on scales orders of magnitude apart. The
  proposal is uniform or Cauchy, the second being the visiting distribution
  of fast simulated annealing (Szu and Hartley 1987), the `q = 2` member of
  the Tsallis family; the general Tsallis visiting distribution is not
  implemented. The initial temperature is calibrated from the variation of
  the objective unless the caller sets it, so the same problem scaled by a
  million is solved as well. The run returns the best point seen and not
  the last, and `converged` is not inferred from the end of the schedule:
  the stationarity measure is Corana's termination rule, read by
  `crit_stationary()`. `sa_run_r()` is an R twin of the compiled loop; the
  two draw from the generator of R in the same order, so from one seed they
  give the same run and the test needs no tolerance. The compiled loop is
  1.74 times faster on an objective costing 0.7 microseconds and 1.34 times
  on one costing 3.8; on an objective of half a millisecond, the cost of an
  inner fit of a modeling layer, the loop overhead is a quarter of a per
  cent.

* `chain()`, which runs optimizers one after another, each starting where
  the previous one finished: `chain(sa(), lbfgs())` explores and then
  descends, and neither method needs code for the other. It is the second
  wrapper after `multistart()`, and the two combine. Each stage carries its
  own criterion and budgets; the point, the value and `converged` are those
  of the last stage, the work is summed, and the trace has a `stage`
  column.

* The ninth check of `check_optimizer()` runs its `maximize()` comparison
  from the same random stream as the `minimize()` run that it is compared
  with, as the eighth check already did. Without it the check required the
  optimizer to converge tightly enough that two independent runs agree to
  1e-5, which the mesh-shrinking methods do, so it passed for the wrong
  reason on them and failed on a global search.

# optimizers7 0.1.0

## Methods

* `prox_grad()`, the proximal gradient method: an accelerated iteration
  with backtracking and adaptive restart. Its loop is written in R, every
  iteration calling a proximal operator supplied by the caller. On one
  seeded lasso draw acceleration mattered on ill-conditioned problems only:
  at a condition number of 3 the plain iteration was slightly faster, at 55
  it took 4153 iterations against 126, and at 480 it did not converge in
  50000 while the accelerated one took 334.

* `optimizer_bounded()` declares whether a method accepts box bounds beside
  the objective. It is `FALSE` for `prox_grad()`, which takes its constraint
  inside the operator; `check_optimizer()` reads it before testing bounds.

* `nonmonotone()`, a third `line_search` subclass implementing the
  condition of Grippo, Lampariello and Lucidi, and the default for `bb()`.
  With `memory = 0` it is `armijo()` value for value. With the default memory
  on Rosenbrock, `bb()` took 68 iterations and 77 evaluations against 82 and
  186 under `armijo()`, accepting eleven uphill steps out of sixty-seven,
  where Armijo accepts none.

* `bb()` takes its fallback step from `1/max|g|` when the secant pair has no
  usable curvature. Freezing the previous step length traps a short step
  (873 of 945 iterations, on bb2), and restarting at `alpha0` does the same
  wherever `alpha0` is too short (1395 rejected steps in 1521 iterations, on
  a boxed quadratic). `alpha_max`, which the SPG literature uses, gives the
  search a direction of length 1e10, and the run stopped a unit from the
  solution reporting convergence. The curvature test is relative,
  `s'y > curv_tol*||s||*||y||`, as in `bfgs()`.

* `start_zeros()` and `start_runif()`, so that a starting value need not be
  written out, and `multistart()` runs its starts in parallel by default.

## Stopping and reporting

* The default gradient tolerance is `1e-6`, in `crit_grad()` and in every
  method that uses it. A line search accepts a step only when the objective
  decreases by a definite amount, so the smallest gradient that a run can
  reach is about `sqrt(2*lambda*eps*|f*|)`, which grows with the value at
  the solution. Every problem in `test_problems()` has `f* = 0`, where that
  floor is zero; a log-likelihood divided by the number of observations is
  of order one at its maximum, and the floor is then around 1e-8. Across
  families, methods and seeds it reached 1.06e-8.

* `gd()`, `cg()`, `bb()`, `newton()`, `bfgs()` and `lbfgs()` no longer
  restate the tolerances in their defaults, so a change to `crit_grad()`
  reaches them.

* A run evaluates its stopping rule before reporting a line-search failure.
  A search that found no acceptable step ended the loop with
  `converged = FALSE` even when the start was already the answer, as it is
  for a caller who passes a closed-form estimate. The rule is evaluated
  without a previous value, so a rule on a change of the objective returns
  `FALSE` and only the state at the point can end the run.

* The gradient check is skipped at a stationary start. Differencing the
  objective along a direction without slope gives the truncation error of
  the difference, so the check warned exactly when the caller supplied the
  optimum.

* A supplied gradient is checked once against the objective, with one
  central difference along the gradient direction at `par`, and a large
  disagreement gives a warning that names both rates. Disable it with
  `options(optimizers7.check_gradient = FALSE)`.

* `newton()` without `gr` and `he` warns when the evaluation budget admits
  fewer than two iterations, and `max_eval` defaults to `Inf`.

## Interface

* Box constraints are two vectors, `lower` and `upper`, instead of a list
  of pairs, so `lower = 0` makes every parameter positive.

* `adam(resample = )` and `finite_sum()` are removed. An objective that
  draws its own minibatches is a closure written in one line, and the two
  required a second objective class, a rule for the criteria that they
  permitted, a name in the criterion machinery, a branch in the compiled
  loop and an argument to select the path.

* `cpp_objective()` is removed. The pointer type
  `double(*)(const arma::vec&)` has no closure, so a statistical objective
  needed C++ globals, which is not reentrant.

* `multistart()` runs its starts sequentially under `pkgload`. Its PSOCK
  workers load the installed copy of the package, so an S7 optimizer built
  in a development namespace dispatched against the methods of the installed
  namespace and the inner runs returned wrong results without an error.

## Validation and documentation

* `check_optimizer()` separates what an optimizer reports from its
  performance: twelve checks on what it reports (that `value` is the
  objective at `par`, that `converged` follows the rule and is not inferred
  from the end of the run, that bounds hold strictly, and others), followed
  by a battery of standard problems whose gaps are printed as information.
  It found two defects in `bundle()`, one of them a stationarity measure
  that contained the trust parameter, which the algorithm shrinks for
  reasons of its own.

* `check_bounds()`, `check_criterion()` and `as_objective()` are exported,
  because a user-defined optimizer cannot honor box constraints or reject a
  stopping rule that it cannot evaluate without them.

* Every constructor cites the paper of its method, every exported topic
  carries a `\value` and an executable example, and the class pages live at
  `X-class` with the plain name as an alias, since two topics that differ
  only in case are one file on Windows.

* A vignette on extending the package, a README with badges, a pkgdown site
  and continuous integration on five platforms.
