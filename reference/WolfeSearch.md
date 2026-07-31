# S7 Class for the Strong Wolfe Line Search

The class
[`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
instantiates.

## Usage

``` r
WolfeSearch(
  label = character(0),
  c1 = integer(0),
  c2 = integer(0),
  max_step = integer(0)
)
```

## Arguments

- c1:

  The sufficient-decrease constant.

- c2:

  The curvature constant.

- max_step:

  The most trial steps allowed.

## Value

An S7 object inheriting from
[`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md).

## See also

[`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
