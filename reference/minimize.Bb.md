# Minimize by Barzilai-Borwein

Runs [`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md)
on the objective: the step length is formed from the previous secant
pair and offered to the line search unaltered, which backtracks only if
it is rejected. It shares the descent loop with
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) and
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md).

## Arguments

- optimizer:

  A `Bb` object.

- fn, par, gr, hess, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `hess` is accepted and ignored. Bounds are taken and removed by
  reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose trace, when kept, reports `bb curvature reset` where a secant pair
had no usable curvature, `bb step clamped` where the step was clamped to
`[alpha_min, alpha_max]`, and `step shortened` where the line search
backtracked.
