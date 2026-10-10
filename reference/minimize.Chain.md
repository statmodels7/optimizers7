# Minimize by a Sequence of Optimizers

Runs the stages in order, passing to each the point that the previous
one reached, and assembles one result from them. The objective, its
derivatives and the bounds are passed to every stage unchanged; only the
starting point changes.

## Arguments

- optimizer:

  A `Chain` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  Passed to every stage as given, with `par` replaced by the answer of
  the previous stage.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
with the point, value, gradient, iteration count and verdict of the last
stage, the summed evaluation counts, the stacked trace and the elapsed
time of the whole chain. `seed` is that of the **first** stage, since a
repeat of the chain has to start from that state.

## Details

Two things happen before the first stage runs. The stopping rule of
every stage is checked by
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md),
so a chain whose last stage cannot evaluate its own rule fails without
running the earlier ones. And the gradient check is switched off for the
stages, since the generic has already made it at the `par` of the
caller; without that, the same warning would be printed once per stage.
