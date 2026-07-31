# S7 Class for Multi-Start

The class
[`multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
instantiates.

## Usage

``` r
MultiStart(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  optimizer = NULL,
  n = integer(0),
  starts = NULL,
  spread = integer(0),
  ncores = NULL,
  distinct_tol = integer(0)
)
```

## Arguments

- optimizer:

  The inner optimiser, run from each starting point.

- n:

  How many starts.

- starts:

  An optional matrix of starting points.

- spread:

  How widely the random starts are scattered.

- ncores:

  How many processes the starts are spread over.

- distinct_tol:

  Objective values closer than this count as one optimum.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
