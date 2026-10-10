# Minimize by Simulated Annealing

Runs [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md)
on the objective: a Metropolis walk over one coordinate at a time, with
the step of each coordinate adapted to hold its acceptance rate near the
target, and the temperature cooled geometrically. Records the state of
the random number generator before the first draw, so the run can be
repeated.

## Arguments

- optimizer:

  An `Sa` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `gr` and `he` are accepted and ignored, since the method uses no
  derivative; rejecting them would force calling code to branch on the
  algorithm. Bounds are taken and removed by reparametrization, as in
  the other methods except
  [`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md).

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
holding the best point seen, with `seed` filled in and `converged`
`TRUE` only when the stopping rule fired (by default, Corana's rule
through
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)).
