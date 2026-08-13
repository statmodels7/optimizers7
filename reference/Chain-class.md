# S7 Class for a Sequence of Optimizers

The class
[`chain`](https://statmodels7.github.io/optimizers7/reference/chain.md)
instantiates.

## Usage

``` r
Chain(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  stages = list()
)
```

## Arguments

- stages:

  The optimizers, in the order they run.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`chain`](https://statmodels7.github.io/optimizers7/reference/chain.md)
