# S7 Class for Armijo Backtracking

The class
[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
instantiates.

## Usage

``` r
ArmijoSearch(
  label = character(0),
  c1 = integer(0),
  shrink = integer(0),
  max_step = integer(0)
)
```

## Arguments

- c1:

  The sufficient-decrease constant.

- shrink:

  The factor the step is multiplied by on each backtrack.

- max_step:

  The most backtracks allowed.

## Value

An S7 object inheriting from
[`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md).

## See also

[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
