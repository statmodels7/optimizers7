# Minimise by the Proximal Bundle Method

Runs
[`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
on the objective.

## Arguments

- optimizer:

  A `Bundle` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `gr` should return a subgradient; `he` is ignored.

## Value

An
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md),
whose `gradient` is the *aggregate* subgradient, the one quantity of
that shape which goes to zero at a solution.
