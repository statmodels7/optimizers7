# Minimize From Many Starting Points

Runs
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
on the objective: generate or take the starting points, run the inner
optimizer from each, and assemble the best.

## Arguments

- optimizer:

  A `MultiStart` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `par` is the first starting point and is used as given; the rest are
  generated around it.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md):
the best run, carrying its own `par`, `value`, `gradient` and
`converged`, with the per-start summary in `trace` and the counts in
`message`. `seed` is the state the whole run began from.

## Details

Three things happen before the first start runs. The objective and the
bounds are validated **once**, so a bad criterion or a start outside its
box is refused before \\n\\ runs are launched rather than \\n\\ times
over. The generator state is captured for the result's `seed`. And the
gradient-consistency check is switched off for the duration, the generic
having already made it at the caller's `par`; without that the same
warning would print once per start.
