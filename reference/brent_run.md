# The Two Phases of Brent's Method

This is the body of
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
for a
[Brent](https://statmodels7.github.io/optimizers7/reference/Brent-class.md)
optimizer: the bracketing phase and the reduction phase, on the
unconstrained scale where bounds are given.

## Usage

``` r
brent_run(optimizer, spec, par, bounds)
```

## Arguments

- optimizer:

  A `Brent` object.

- spec:

  The objective handle from
  [`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md).

- par:

  The starting value, a single number.

- bounds:

  What
  [`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
  returned.

## Value

A list in the shape that
[`build_result()`](https://statmodels7.github.io/optimizers7/reference/build_result.md)
consumes.
