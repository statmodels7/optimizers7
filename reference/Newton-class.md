# S7 Class for Newton's Method

The class
[`newton`](https://statmodels7.github.io/optimizers7/reference/newton.md)
instantiates.

## Usage

``` r
Newton(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  line_search = NULL,
  hessian_mod = character(0),
  floor = integer(0)
)
```

## Arguments

- step, line_search:

  As in
  [`newton`](https://statmodels7.github.io/optimizers7/reference/newton.md).

- hessian_mod:

  How an indefinite Hessian is repaired.

- floor:

  The smallest eigenvalue the repaired Hessian may have.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`newton`](https://statmodels7.github.io/optimizers7/reference/newton.md)
