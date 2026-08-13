# S7 Class for Simulated Annealing

The class
[`sa`](https://statmodels7.github.io/optimizers7/reference/sa.md)
instantiates.

## Usage

``` r
Sa(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  visiting = character(0),
  t0 = NULL,
  cooling = integer(0),
  cycles = integer(0),
  steps = integer(0),
  step = integer(0),
  target_accept = integer(0),
  adjust = integer(0),
  n_eps = integer(0)
)
```

## Arguments

- visiting:

  Which proposal the walk uses.

- t0:

  The initial temperature, or `NULL` to calibrate it.

- cooling:

  The geometric cooling factor.

- cycles, steps:

  The work done at each temperature.

- step:

  The initial step, relative to the starting value.

- target_accept:

  The acceptance rate the step adaptation aims at.

- adjust:

  How hard the step is adjusted towards that rate.

- n_eps:

  How many temperature levels the stopping rule looks back over.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`sa`](https://statmodels7.github.io/optimizers7/reference/sa.md)
