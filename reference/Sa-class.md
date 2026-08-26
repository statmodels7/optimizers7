# S7 Class for Simulated Annealing

An optimizer holding the annealing schedule and the step adaptation:
which proposal the walk draws from, the temperature and how it cools,
how much work is done at each level, and the acceptance rate the
per-coordinate steps are tuned to. Built by
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md). It
is the package's only global search and the only method whose result
depends on the state of the random number generator.

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

  `"uniform"` or `"cauchy"`, which proposal the walk draws its moves
  from.

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

  The acceptance rate the step adaptation aims at.

- adjust:

  How hard the step is adjusted towards that rate.

- n_eps:

  How many temperature levels the stopping rule looks back over.

## Value

An S7 object of class `Sa` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the nine properties above beside the seven shared ones.

## Details

Beyond the seven properties every optimizer has, an `Sa` carries nine of
its own. They fall into three groups: the proposal (`visiting`, `step`),
the schedule (`t0`, `cooling`, `cycles`, `steps`), and the adaptation
and its stopping rule (`target_accept`, `adjust`, `n_eps`).

## See also

[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) for
the constructor,
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
for handing its answer to a local method.
