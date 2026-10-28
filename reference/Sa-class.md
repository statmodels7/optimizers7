# S7 Class for Simulated Annealing

This class is an optimizer that holds the annealing schedule and the
step adaptation: the proposal from which the walk draws, the temperature
and how it cools, how much work is done at each level, and the
acceptance rate to which the per-coordinate steps are tuned. It is built
by [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md).
Its result depends on the state of the random number generator, which
the result records in `seed`.

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
  proposal = character(0),
  t0 = NULL,
  cooling = integer(0),
  cycles = integer(0),
  steps = integer(0),
  step = integer(0),
  target_accept = integer(0),
  adjust = integer(0),
  stationarity_levels = integer(0)
)
```

## Arguments

- proposal:

  `"uniform"` or `"cauchy"`, the proposal from which the walk draws its
  moves.

- t0:

  The initial temperature, or `NULL` to calibrate it from the objective.

- cooling:

  The geometric cooling factor, in \\(0, 1)\\.

- cycles, steps:

  The work done at each temperature level: `cycles * steps * p`
  evaluations.

- step:

  The initial step, relative to the starting value.

- target_accept:

  The target acceptance rate of the step adaptation.

- adjust:

  How strongly the step is adjusted towards that rate.

- stationarity_levels:

  How many temperature levels the stationarity measure covers.

## Value

An S7 object of class `Sa` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the nine properties above beside the seven shared ones.

## Details

Beyond the seven properties of every optimizer, an `Sa` carries nine of
its own. They fall into three groups: the proposal (`proposal`, `step`),
the schedule (`t0`, `cooling`, `cycles`, `steps`), and the adaptation
and its stopping rule (`target_accept`, `adjust`,
`stationarity_levels`).

## See also

[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) for
the constructor,
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
for handing its answer to a local method.
