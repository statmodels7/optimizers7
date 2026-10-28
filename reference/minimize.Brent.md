# Minimize by Brent's Method

Runs
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
on a function of one variable: bracket the minimum by golden expansion,
then reduce the bracket by parabolic and golden-section steps until the
stopping rule fires or the bracket reaches the resolution of the
iterate.

## Arguments

- optimizer:

  A `Brent` object.

- fn, par, gr, hess, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `gr` and `hess` are accepted and ignored, since the method uses no
  derivative. Bounds are taken and removed by reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose `gradient` is `NULL` and whose trace has the columns `iteration`,
`value`, `stationarity` (the half-width of the bracket) and `safeguard`,
which names the step taken at each iteration: `parabolic` or `golden`.
