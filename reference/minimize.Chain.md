# Minimize by a Sequence of Optimizers

Runs the stages in order, handing each the point the previous one
reached, and assembles one result from them. The objective, its
derivatives and the bounds are passed to every stage unchanged; only the
starting point moves.

## Arguments

- optimizer:

  A `Chain` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  Passed to every stage as given, with `par` replaced by the previous
  stage's answer.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
carrying the last stage's point, value, gradient, iteration count and
verdict, the summed evaluation counts, the stacked trace and the elapsed
time of the whole chain. `seed` is the **first** stage's, that being the
state a repeat of the chain has to start from.

## Details

Two things happen before the first stage runs. Every stage's stopping
rule is put to
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md),
so a chain whose last stage cannot evaluate its own rule fails without
spending the earlier ones. And the gradient-consistency check is
switched off for the duration, the generic having already made it once
at the caller's `par`; without that, the same warning would print once
per stage.
