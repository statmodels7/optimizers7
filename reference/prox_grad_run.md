# Run the Proximal Gradient Loop

The iteration behind
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md):
a backtracked gradient step on the smooth part, the proximal operator
applied to its result, and the momentum extrapolation with its restart.
It counts its own evaluations, including those of a gradient obtained by
differences.

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

A list in the shape that
[`build_result()`](https://statmodels7.github.io/optimizers7/reference/build_result.md)
consumes: the point, the total objective there, the mapping, the counts,
the iteration count, the verdict and the trace.

## Details

The loop is written in R, as that of
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md).
Every iteration calls three R functions supplied by the caller: the
objective, its gradient and the proximal operator.

The proximal gradient mapping is read at the iterate, the point that is
reported, and not at the extrapolated point. The trace records
`step shortened` at an iteration whose step was backtracked and
`momentum restart` at one whose momentum was reset.
