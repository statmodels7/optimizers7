# Package index

## Running an optimization

One generic for every algorithm. The kind of objective is handled
separately, by
[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md),
so an algorithm is written once.

- [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
  : Minimize a Function
- [`maximize()`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
  : Maximize a Function
- [`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
  : S7 Class for the Result of an Optimization

## Smooth problems

Methods that use derivatives of the objective. They are fast when the
objective is smooth, and they cannot certify a solution at a kink.

- [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
  : Newton's Method with a Modified Hessian
- [`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
  : Chebyshev's Third-Order Method
- [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  : BFGS
- [`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
  : Limited-Memory BFGS
- [`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md) :
  Nonlinear Conjugate Gradients
- [`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) :
  The Barzilai-Borwein Method
- [`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) :
  Gradient Descent

## Noisy and large problems

Adam, which tolerates a gradient that points downhill only on average.
It draws no minibatches of its own: an objective that resamples is a
closure written by the caller, who defines what an observation is.

- [`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
  : Adaptive Moment Estimation

## Problems with a kink

Methods for objectives with kinks, at which a descent method reaches the
solution without being able to certify it, because the subgradient that
it evaluates does not become small.
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
tests whether zero lies in the convex hull of its subgradients,
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
handles the kink through a proximal operator, and
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
and
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
use no derivative.

- [`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
  : The Proximal Bundle Method
- [`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
  : Proximal Gradient Method
- [`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
  : The Nelder-Mead Simplex Method
- [`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
  : Pattern Search, With Coordinate or Random Polling

## One parameter

A function of one variable, minimized without derivatives by bracketing
the minimum and reducing the bracket.

- [`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
  : Brent's Method in One Dimension

## Searching before descending

A global search for problems on which the answer depends on the basin in
which the run begins. Its result is meant to be passed to a local
method, which finds the minimum.

- [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) :
  Simulated Annealing With an Adaptive Step

## Wrapping an optimizer

Both are themselves optimizers, so they combine with every method above
and with each other. A multi-start reports the number of distinct
optimal values reached from its starts; a chain runs its stages in
order, and each stage starts where the previous one finished.

- [`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
  : Run an Optimizer From Many Starting Points
- [`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
  : Run Optimizers One After Another

## Stopping rules

A criterion is an object that implements one generic, so a user-defined
rule is treated like a shipped one. An optimizer rejects a rule that it
cannot evaluate instead of accepting one that never fires.

- [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  : S7 Class for Convergence Criteria
- [`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
  : Has the Stopping Rule Been Met?
- [`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)
  : What a Criterion Needs From the Iteration
- [`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
  : Stop When the Gradient Is Small
- [`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)
  : Stop When the Objective Stops Moving (Absolute)
- [`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
  : Stop When the Objective Stops Moving (Relative)
- [`crit_abs_par()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md)
  : Stop When the Parameters Stop Moving (Absolute)
- [`crit_rel_par()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_par.md)
  : Stop When the Parameters Stop Moving (Relative)
- [`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
  : Stop When the Method's Own Measure of Progress Is Small
- [`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md)
  : The Rule That Never Fires
- [`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
  : Stop When Any of Several Rules Fires
- [`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)
  : Stop Only When Every Rule Fires

## Line searches

How far to go along a direction, once the direction is chosen. The
nonmonotone search is the default of the Barzilai-Borwein method, whose
step length is the method itself and which relies on steps that increase
the objective.

- [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  : S7 Class for Line Searches
- [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
  : Backtracking Line Search with the Armijo Condition
- [`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
  : Line Search Satisfying the Strong Wolfe Conditions
- [`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
  : Nonmonotone Backtracking

## Starting values

A starter takes the place of the vector of starting values and is
resolved on the unconstrained scale, so a single constant gives a
sensible value for every kind of parameter and every draw lies inside
its bounds.

- [`starter()`](https://statmodels7.github.io/optimizers7/reference/starter-class.md)
  : S7 Class for a Starting-Value Generator
- [`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
  : Start From Zero on the Unconstrained Scale
- [`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md)
  : Start From a Uniform Draw on the Unconstrained Scale
- [`starting_values()`](https://statmodels7.github.io/optimizers7/reference/starting_values.md)
  : Produce a Vector of Starting Values
- [`infer_npar()`](https://statmodels7.github.io/optimizers7/reference/infer_npar.md)
  : How Many Parameters the Objective Takes

## The objective

A function of the parameter vector. The generic is the extension point:
a caller with another kind of objective registers one method and every
algorithm accepts that objective.

- [`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md)
  : Normalize an Objective for the Optimizers

## Writing a new optimizer

The requirements on a new
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
method, and the functions that help to meet them.
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
runs twelve checks on what an optimizer reports, and reports separately
the gaps that it reaches on the test problems, and these gaps do not
enter the verdicts.

- [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  : S7 Class for Optimization Algorithms
- [`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
  : Reject a Stopping Rule That an Optimizer Cannot Evaluate
- [`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
  : Normalize Box Constraints
- [`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
  : What an Optimizer Can Offer a Stopping Rule
- [`optimizer_bounded()`](https://statmodels7.github.io/optimizers7/reference/optimizer_bounded.md)
  : Whether an Optimizer Takes Box Bounds
- [`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)
  : The Bound Transform, and Its First Two Derivatives
- [`bounded_forward()`](https://statmodels7.github.io/optimizers7/reference/bounded_forward.md)
  : The Forward Bound Transform
- [`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
  : Check What an Optimizer Reports
- [`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
  : The Standard Test Problems
