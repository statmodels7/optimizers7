# S7 Class for Adam

The class
[`adam`](https://statmodels7.github.io/optimizers7/reference/adam.md)
instantiates.

## Usage

``` r
Adam(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  alpha = integer(0),
  beta1 = integer(0),
  beta2 = integer(0),
  eps = integer(0),
  decay = integer(0),
  amsgrad = logical(0)
)
```

## Arguments

- alpha:

  The learning rate.

- beta1, beta2:

  Decay rates for the two moment estimates.

- eps:

  The denominator floor.

- decay:

  Rate at which the learning rate is reduced.

- amsgrad:

  Whether to hold the second moment at its running maximum.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`adam`](https://statmodels7.github.io/optimizers7/reference/adam.md)
