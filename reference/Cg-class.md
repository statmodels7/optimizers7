# S7 Class for Conjugate Gradients

The class
[`cg`](https://statmodels7.github.io/optimizers7/reference/cg.md)
instantiates.

## Usage

``` r
Cg(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  beta = character(0),
  restart_every = integer(0),
  step = integer(0),
  line_search = NULL
)
```

## Arguments

- beta:

  Which update formula.

- restart_every:

  How often the method is restarted at steepest descent.

- step:

  The initial step length offered to the line search.

- line_search:

  A
  [`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`cg`](https://statmodels7.github.io/optimizers7/reference/cg.md)
