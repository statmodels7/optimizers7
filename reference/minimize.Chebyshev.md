# Minimize by Chebyshev's Method

Minimize by Chebyshev's Method

## Arguments

- optimizer:

  A `Chebyshev` object.

- fn, par, gr, he, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

- t3:

  `NULL`, or a function of the parameters and a direction returning the
  third derivative of `fn` contracted twice with the direction, a vector
  as long as the parameters.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).
