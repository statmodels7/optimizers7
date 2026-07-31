# What a Derivative-Free Method Can Offer a Stopping Rule

The objective and a stationarity measure, but no gradient.

## Arguments

- optimizer:

  A `NelderMead` or `Compass` object.

## Value

A character vector.

## Details

A rule reading a gradient is refused rather than left testing `NULL` at
every iteration and never firing.
[`crit_stationary`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
is what takes its place.
