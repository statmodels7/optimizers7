# Minimize by Pattern Search

Runs
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
on the objective: poll the directions around the current point at the
current radius, move to an improvement if one is found and expand, and
shrink the radius when the poll fails everywhere.

## Arguments

- optimizer:

  A `Compass` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `gr` and `he` are accepted and ignored, the method using no
  derivative. Bounds are taken and removed by reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose `gradient` is `NULL` and whose trace carries a `stationarity`
column holding the poll size. With `directions = "mads"` the `seed` is
recorded, the poll drawing from R's generator; with `"coordinate"` it is
`NULL` and the run repeats without one.
