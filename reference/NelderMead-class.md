# S7 Class for Nelder-Mead

The class
[`nelder_mead`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
instantiates.

## Usage

``` r
NelderMead(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  adaptive = logical(0),
  max_restarts = integer(0),
  degenerate_tol = integer(0),
  simplex = NULL
)
```

## Arguments

- step:

  Relative size of the initial simplex.

- adaptive:

  Whether to use dimension-dependent coefficients.

- max_restarts:

  How many times a degenerate simplex may be rebuilt.

- degenerate_tol:

  The conditioning below which it is rebuilt.

- simplex:

  An optional starting simplex.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`nelder_mead`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
