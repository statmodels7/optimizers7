# Minimize by Barzilai-Borwein

Runs [`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md)
on the objective: form the step length from the previous secant pair,
offer that step to the line search unaltered, and backtrack only if it
is rejected. Shares the descent loop with
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) and
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md).

## Arguments

- optimizer:

  A `Bb` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `he` is accepted and ignored. Bounds are taken and removed by
  reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose trace, when kept, reports `bb curvature reset` where a secant pair
carried none and `step shortened` where the step was clamped or
backtracked.
