# S7 Class for Pattern Search

The class
[`compass`](https://statmodels7.github.io/optimizers7/reference/compass.md)
instantiates.

## Usage

``` r
Compass(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  directions = character(0),
  opportunistic = logical(0),
  expand = integer(0),
  shrink = integer(0)
)
```

## Arguments

- step:

  Initial poll size, relative to the starting value.

- directions:

  Either `"mads"` or `"coordinate"`.

- opportunistic:

  Whether to accept the first improvement found.

- expand, shrink:

  Factors applied to the poll size.

## Value

An S7 object inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## See also

[`compass`](https://statmodels7.github.io/optimizers7/reference/compass.md)
