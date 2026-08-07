# Package index

## Running an optimization

One generic for every algorithm. Which of the three shapes the objective
arrived in is settled separately, so an algorithm is written once.

- [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
  : Minimize a Function
- [`maximize()`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
  : Maximize a Function
- [`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
  : S7 Class for the Result of an Optimization

## Smooth problems

Methods that build a model of the surface from derivatives. Fast when
the objective is smooth, and misled when it is not.

- [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
  : Newton's Method with a Modified Hessian
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
closure, written by the caller, who knows what an observation is.

- [`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
  : Adaptive Moment Estimation

## Problems with a kink

Where a descent method arrives at the answer and cannot certify it,
because the subgradient it evaluates never becomes small.

- [`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
  : The Proximal Bundle Method
- [`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
  : The Nelder-Mead Simplex Method
- [`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
  : Pattern Search, With Coordinate or Random Polling

## Wrapping an optimizer

Multi-start is itself an optimizer, so it composes with everything
above. What it reports that a single run cannot is the number of
distinct optima.

- [`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
  : Run an Optimizer From Many Starting Points

## Stopping rules

A criterion is an object implementing one generic, so a user-defined
rule is treated like a shipped one. An optimizer rejects a rule it
cannot evaluate rather than accepting one that never fires.

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
  : Never Stop Early
- [`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
  : Stop When Any of Several Rules Fires
- [`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)
  : Stop Only When Every Rule Fires

## Line searches

How far to go along a direction, once the direction is chosen. The
non-monotone rule is not an exotic option: a method whose step length is
itself the method, as Barzilai-Borwein’s is, needs a search that does
not reshape it.

- [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  : S7 Class for Line Searches
- [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
  : Backtracking Line Search with the Armijo Condition
- [`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
  : Line Search Satisfying the Strong Wolfe Conditions
- [`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
  : Nonmonotone Backtracking

## Starting values

A starter stands in for the vector of starting values, and is resolved
on the unconstrained scale, so a single constant means something
sensible for every kind of parameter and no draw can land outside its
bounds.

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
a caller with its own kind of objective registers one method and every
algorithm accepts it.

- [`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md)
  : Normalize an Objective for the Optimizers

## Writing an optimizer of your own

What a
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
method of your own has to promise, and the pieces that let it.
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
tests the contract; the battery under it reports the power; they are
separate because they are different questions.

- [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  : S7 Class for Optimization Algorithms
- [`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
  : Refuse a Stopping Rule an Optimizer Cannot Evaluate
- [`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
  : Normalize Box Constraints
- [`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
  : What an Optimizer Can Offer a Stopping Rule
- [`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)
  : The Bound Transform, and Its First Two Derivatives
- [`bounded_forward()`](https://statmodels7.github.io/optimizers7/reference/bounded_forward.md)
  : The Forward Bound Transform
- [`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
  : Check That an Optimizer Keeps Its Promises
- [`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
  : The Standard Test Problems
