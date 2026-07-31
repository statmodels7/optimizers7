# S7 Class for the Barzilai-Borwein Method

The class
[`bb`](https://statmodels7.github.io/optimizers7/reference/bb.md)
instantiates.

## Usage

``` r
Bb(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  variant = character(0),
  alpha0 = integer(0),
  alpha_min = integer(0),
  alpha_max = integer(0),
  step = integer(0),
  line_search = NULL
)
```

## Arguments

- variant:

  Which step-length formula.

- alpha0:

  The step length used before there is a secant pair.

- alpha_min, alpha_max:

  Bounds on it.

- step:

  The initial multiplier offered to the line search.

- line_search:

  A
  [`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`bb`](https://statmodels7.github.io/optimizers7/reference/bb.md)
