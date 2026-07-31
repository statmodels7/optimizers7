# S7 Class for Limited-Memory BFGS

The class
[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
instantiates.

## Usage

``` r
Lbfgs(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  line_search = NULL,
  memory = integer(0),
  curv_tol = integer(0)
)
```

## Arguments

- step, line_search:

  As in
  [`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md).

- memory:

  How many secant pairs to keep.

- curv_tol:

  The curvature threshold below which a pair is not stored.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
