# Minimise by Pattern Search

Runs
[`compass`](https://statmodels7.github.io/optimizers7/reference/Compass.md)
on the objective.

## Arguments

- optimizer:

  A `Compass` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `gr` and `he` are accepted and ignored.

## Value

An
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).
