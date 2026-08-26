# Minimize by Conjugate Gradients

Runs [`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md)
on the objective: bend each direction with the previous one by the
chosen \\\beta\\, and let the line search choose how far along it to go.
Shares the descent loop with
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) and
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md).

## Arguments

- optimizer:

  A `Cg` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `he` is accepted and ignored. Bounds are taken and removed by
  reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose trace, when kept, reports `cg restart` at each iteration where the
bend was clamped to zero or the direction was replaced by the gradient.
