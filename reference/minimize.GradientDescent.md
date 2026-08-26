# Minimize by Gradient Descent

Runs [`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md)
on the objective: take \\-g\\ as the direction and let the line search
choose how far. Shares the descent loop with
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md) and
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md),
which differ from it only in how the direction is formed.

## Arguments

- optimizer:

  A `GradientDescent` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `he` is accepted and ignored. Bounds are taken and removed by
  reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose `gradient` is \\\nabla f\\ at `par` and whose trace, when kept,
carries a `gnorm` column.
