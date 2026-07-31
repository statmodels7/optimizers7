# An Ordinary R Function as an Objective

The common case: `fn(par)` returns a number and, if supplied, `gr(par)`
returns the gradient.

## Arguments

- fn:

  A function of the parameter vector.

- gr:

  An optional gradient function, or `NULL` for finite differences.

- ...:

  Unused.

## Value

An objective handle; see
[`as_objective`](https://statmodels7.github.io/optimizers7/reference/as_objective.md).
