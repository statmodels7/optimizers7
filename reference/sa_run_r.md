# The Annealing Loop in R

The loop `sa_run()` replaces, kept so the compiled route has something
to be compared against that shares none of its code.

## Usage

``` r
sa_run_r(
  fn,
  par,
  cauchy = FALSE,
  t0 = -1,
  cooling = 0.85,
  cycles = 3,
  steps = 10,
  step = 1,
  target_accept = 0.5,
  adjust = 2,
  n_eps = 4,
  maxit = 100
)
```

## Arguments

- fn:

  The objective.

- par:

  The starting value.

- cauchy:

  Whether the proposal is Cauchy rather than uniform.

- t0:

  The initial temperature, or a non-positive value to calibrate it.

- cooling, cycles, steps, step, target_accept, adjust, n_eps, maxit:

  As in
  [`sa`](https://statmodels7.github.io/optimizers7/reference/sa.md).

## Value

A list with `par`, `value` and `n_value`.

## Details

It draws from R's generator in the same order as the kernel, so from one
seed the two runs are the same run and the comparison needs no
tolerance. The order is what has to match: one uniform per proposal, and
a second one for the Metropolis test ONLY when the proposal is uphill,
which is where a transcription most easily drifts.
