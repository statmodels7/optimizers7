# S7 Class for the Proximal Gradient Method

The class
[`prox_grad`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
instantiates.

## Usage

``` r
ProxGrad(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  prox = function() NULL,
  g = function() NULL,
  accelerate = logical(0),
  step = integer(0),
  shrink = integer(0),
  restart = logical(0)
)
```

## Arguments

- prox:

  The proximal operator of the non-smooth part.

- g:

  The value of the non-smooth part.

- accelerate:

  Whether the momentum extrapolation is applied.

- step:

  The initial step length.

- shrink:

  The backtracking factor.

- restart:

  Whether an increase in the objective resets the momentum.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`prox_grad`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
