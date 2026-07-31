# Minimise From Many Starting Points

Runs
[`multistart`](https://statmodels7.github.io/optimizers7/reference/MultiStart.md)
on the objective.

## Arguments

- optimizer:

  A `MultiStart` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Value

An
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md):
the best run, with the per-start summary in its `trace`.
