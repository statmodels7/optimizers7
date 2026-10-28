# Minimize by Gradient Descent

Runs [`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md)
on the objective: the direction is \\-g\\, and the line search chooses
the step length. It shares the descent loop with
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md) and
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md),
which differ from it only in the direction.

## Arguments

- optimizer:

  A `GradientDescent` object.

- fn, par, gr, hess, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `hess` is accepted and ignored. Bounds are taken and removed by
  reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose `gradient` is \\\nabla f\\ at `par` and whose trace, when kept,
carries a `gnorm` column.
