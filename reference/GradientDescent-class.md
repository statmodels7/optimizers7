# S7 Class for Gradient Descent

An optimizer holding the initial step length and the line search, and
nothing else: steepest descent carries no model of the surface, so there
is nothing else to hold. Built by
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md).

## Usage

``` r
GradientDescent(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  line_search = NULL
)
```

## Arguments

- step:

  The initial step length offered to the line search.

- line_search:

  A
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object.

## Value

An S7 object of class `GradientDescent` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with `step` and `line_search` beside the seven shared properties.

## Details

Beyond the seven properties every optimizer has, a `GradientDescent`
carries `step` and `line_search`, both shared with the other line-search
methods. It is the smallest optimizer class in the package.

## See also

[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) for
the constructor,
[Cg](https://statmodels7.github.io/optimizers7/reference/Cg-class.md)
and
[Bb](https://statmodels7.github.io/optimizers7/reference/Bb-class.md)
for the two first-order methods that carry a little more.
