# The Annealing Loop in R

The annealing loop written in R, the twin of the compiled `sa_run()`. It
exists so that the kernel has something to be compared against that
shares none of its code, and it is used by the tests alone.

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

  The objective, a function of the parameter vector.

- par:

  The starting value, a numeric vector.

- cauchy:

  `TRUE` for the Cauchy proposal, `FALSE` for the uniform one.

- t0:

  The initial temperature, or a non-positive value to calibrate it. `-1`
  is what
  [`minimize.Sa()`](https://statmodels7.github.io/optimizers7/reference/minimize.Sa.md)
  passes for `t0 = NULL`.

- cooling, cycles, steps, step, target_accept, adjust, n_eps, maxit:

  As in
  [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md).

## Value

A list with `par` (the best point seen), `value` (the objective there)
and `n_value` (the evaluation count). It reports no trace and no
stationarity measure, the comparison being of the walk alone.

## Details

It draws from R's generator in the same order as the kernel, so from one
seed the two are the same run and the comparison needs no tolerance at
all. The order is the thing that has to match: one uniform per proposal,
and a second one for the Metropolis test only when the proposal is
uphill. A transcription that consumes a uniform unconditionally passes
every test written on the answer and fails this one immediately.
