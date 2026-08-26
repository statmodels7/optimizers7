# Run the Proximal Gradient Loop

The iteration behind
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md):
a backtracked gradient step on the smooth part, the proximal operator
applied to its result, and the momentum extrapolation with its restart.
Counts its own evaluations, so the caller need not.

## Usage

``` r
prox_grad_run(optimizer, spec, par)
```

## Arguments

- optimizer:

  A `ProxGrad` object.

- spec:

  The objective handle from
  [`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md).

- par:

  The starting point, a numeric vector.

## Value

A list in the shape
[`build_result()`](https://statmodels7.github.io/optimizers7/reference/build_result.md)
consumes: the point, the total objective there, the mapping, the counts,
the iteration count, the verdict and the trace.

## Details

This is the one method in the package written in R instead of compiled.
Every iteration calls the objective, its gradient and the proximal
operator, all three R functions supplied by the caller, and the loop
around them costs a fraction of a microsecond against those. Compiling
it would move the callbacks and change nothing else.

The stationarity measure is read **at the iterate** and not at the
extrapolated point. With momentum the two differ, and reading the
extrapolated one leaves a mapping that never vanishes.
