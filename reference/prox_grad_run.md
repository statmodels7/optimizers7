# Run the Proximal Gradient Loop

The iteration behind
[`prox_grad`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md):
a backtracked gradient step on the smooth part, the proximal operator
applied to its result, and the momentum extrapolation with its restart.

## Usage

``` r
prox_grad_run(optimizer, spec, par)
```

## Arguments

- optimizer:

  A `ProxGrad` object.

- spec:

  The objective handle from
  [`as_objective`](https://statmodels7.github.io/optimizers7/reference/as_objective.md).

- par:

  The starting point.

## Value

A list in the shape
[`build_result`](https://statmodels7.github.io/optimizers7/reference/build_result.md)
consumes.

## Details

Written in R rather than compiled, because every iteration calls the
objective, its gradient and the proximal operator, all of which are R
functions supplied by the caller; the loop around them costs a fraction
of a microsecond against those.
