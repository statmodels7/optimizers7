# Minimize by a Sequence of Optimizers

Runs each stage from where the previous one finished.

## Arguments

- optimizer:

  A `Chain` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
  passed to every stage.

## Value

An
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).
