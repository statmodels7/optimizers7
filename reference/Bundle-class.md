# S7 Class for the Proximal Bundle Method

The class
[`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
instantiates.

## Usage

``` r
Bundle(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  t0 = integer(0),
  t_min = integer(0),
  t_max = integer(0),
  m_serious = integer(0),
  bundle_size = integer(0),
  qp_iters = integer(0),
  qp_tol = integer(0)
)
```

## Arguments

- t0:

  Initial proximity weight.

- t_min, t_max:

  Bounds on it.

- m_serious:

  Fraction of the predicted decrease a serious step must achieve.

- bundle_size:

  Largest number of linearizations kept.

- qp_iters, qp_tol:

  Effort spent on the subproblem.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
