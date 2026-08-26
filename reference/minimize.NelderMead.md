# Minimize by Nelder-Mead

Runs
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
on the objective: build a simplex of \\p+1\\ vertices from `par` and
`step`, then reflect, expand, contract or shrink it until its diameter
falls below the tolerance.

## Arguments

- optimizer:

  A `NelderMead` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `gr` and `he` are accepted and ignored, the method using no
  derivative; refusing them would force calling code to branch on the
  algorithm. Bounds are taken and removed by reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose `gradient` is `NULL` and whose trace carries a `stationarity`
column holding the simplex diameter. The `safeguard` column names the
operation taken at each iteration: `reflect`, `expand`, `contract in`,
`contract out`, `shrink`, or `restart` when a degenerate simplex was
rebuilt.

## Details

A `simplex` supplied on the optimizer is checked here against the
starting value, since only now is the number of parameters known: it
must have \\p+1\\ rows and \\p\\ columns, and a mismatch reports both
shapes.
