# Minimize by Conjugate Gradients

Runs [`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md)
on the objective: each direction is combined with the previous one
through the chosen \\\beta\\, and the line search chooses the step
length. It shares the descent loop with
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) and
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md).

## Arguments

- optimizer:

  A `Cg` object.

- fn, par, gr, hess, lower, upper, ...:

  As in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
  `hess` is accepted and ignored. Bounds are taken and removed by
  reparametrization.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
whose trace, when kept, reports `cg restart` at each iteration where
\\\beta\\ was set to zero, because the Polak-Ribière value was negative
or because `restart_every` iterations had passed, and `cg not descent`
at each iteration where the direction was replaced by the negative
gradient.
