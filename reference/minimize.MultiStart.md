# Minimize From Many Starting Points

Runs
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
on the objective: it generates or takes the starting points, runs the
inner optimizer from each, and assembles the result from the best run.

## Arguments

- optimizer:

  A `MultiStart` object.

- fn, par, gr, hess, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  Unless `starts` was given, `par` is the first starting point and is
  used as given, and the others are generated around it. When `starts`
  was given, `par` is used only to check the objective and the bounds
  and to fix the number of parameters.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
with the `par`, `value`, `gradient`, `converged` and `criterion_met` of
the best run, the `counts` summed over the successful starts,
`iterations` equal to the number of starts, the summary of the starts in
`trace` and the counts in `message`. `seed` is the state from which the
whole run began.

## Details

Before the first start runs, the objective, the bounds and `par` are
validated once, so a bad criterion or a `par` outside its box is
rejected before the \\n\\ runs are launched; a row of a user-supplied
`starts` matrix is not checked in advance, and a row outside the box
fails in its own run. The state of the generator is recorded for the
`seed` of the result. The gradient check is switched off for the inner
runs, since the generic has already made it at the `par` of the caller;
without that the same warning would be printed once per start.
