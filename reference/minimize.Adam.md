# Minimize by Adam

Runs
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
on the objective: one gradient per iteration, the two exponentially
weighted moments updated from it, and a coordinatewise step taken with
no line search and no test that the objective fell.

## Arguments

- optimizer:

  An `Adam` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `he` is accepted and ignored. Bounds are taken and removed by
  reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).
Under the default
[`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md)
it reports `converged = FALSE` and `criterion_met`
`iteration budget reached`, the run having ended on its budget.
