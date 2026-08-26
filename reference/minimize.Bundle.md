# Minimize by the Proximal Bundle Method

Runs
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
on the objective: a cutting-plane model built from the subgradients
collected so far, a proximal subproblem solved for the step, and an
acceptance test that turns each trial into a serious step or a null one.

## Arguments

- optimizer:

  A `Bundle` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `gr` should return a **subgradient**, not a difference taken across a
  kink; `he` is accepted and ignored. Bounds are taken and removed by
  reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose `gradient` is the **aggregate** subgradient, the one quantity of
that shape which goes to zero at a solution, and whose `message` reports
the serious and null step counts.
