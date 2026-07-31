# S7 Class for the Nonmonotone Line Search

The class
[`nonmonotone`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
instantiates.

## Usage

``` r
NonmonotoneSearch(
  label = character(0),
  c1 = integer(0),
  shrink = integer(0),
  memory = integer(0),
  max_step = integer(0)
)
```

## Arguments

- c1:

  The sufficient-decrease constant.

- shrink:

  The factor the step is multiplied by on each backtrack.

- memory:

  How many earlier values the reference looks back over.

- max_step:

  The most backtracks allowed.

## Value

An S7 object inheriting from
[`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md).

## See also

[`nonmonotone`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
