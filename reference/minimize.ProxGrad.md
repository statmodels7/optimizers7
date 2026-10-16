# Minimize by the Proximal Gradient Method

Runs
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
on the objective: a backtracked gradient step on the smooth part, the
proximal operator applied to the result, and the momentum extrapolation
with its restart. Box bounds are rejected here, with a message that
names `prox` as the place for the constraint.

## Arguments

- optimizer:

  A `ProxGrad` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  A finite `lower` or `upper` raises an error; `he` is accepted and
  ignored.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose `value` is the total objective \\f(x) + g(x)\\, which is why `g`
is required alongside `prox`, and whose `gradient` is the proximal
gradient mapping at `par` rather than \\\nabla f\\.
