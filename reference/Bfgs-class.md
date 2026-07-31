# S7 Class for BFGS

The class
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
instantiates.

## Usage

``` r
Bfgs(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  line_search = NULL,
  curv_tol = integer(0),
  max_skip = integer(0)
)
```

## Arguments

- step, line_search:

  As in
  [`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md).

- curv_tol:

  The curvature threshold below which the update is skipped.

- max_skip:

  Consecutive skips before the approximation is reset.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
