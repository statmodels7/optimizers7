# Minimize by Simulated Annealing

Runs [`sa`](https://statmodels7.github.io/optimizers7/reference/sa.md)
on the objective.

## Arguments

- optimizer:

  An `Sa` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `gr` and `he` are accepted and ignored: the method uses no derivative,
  and rejecting them would force calling code to branch on the
  algorithm.

## Value

An
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).
