# S7 Class for Armijo Backtracking

This class is a line search that starts from the full step and
multiplies it by `shrink` until the objective decreases by enough.
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
builds it. The search evaluates the objective at trial points and never
the gradient, which makes it the cheaper search per iteration.

## Usage

``` r
ArmijoSearch(
  label = character(0),
  c1 = integer(0),
  shrink = integer(0),
  max_step = integer(0),
  resolution = NULL
)
```

## Arguments

- c1:

  The sufficient-decrease constant.

- shrink:

  The factor by which the step is multiplied on each backtrack.

- max_step:

  The most trial points allowed, including the full step. It is a count,
  not a length.

- resolution:

  The smallest difference in the objective that the objective can
  resolve, a number or a function of no arguments.

## Value

An S7 object of class `ArmijoSearch` inheriting from
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md),
with the four properties above beside `label`.

## Details

Beyond the `label` of every line search, an `ArmijoSearch` holds `c1`,
`shrink`, `max_step` and `resolution`. It has no `c2`, since it has no
curvature condition;
[`line_search_spec()`](https://statmodels7.github.io/optimizers7/reference/line_search_spec.md)
fills that field with a value that the compiled code ignores, so that
the three searches are passed to the compiled code as lists of the same
shape.

## See also

[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
for the constructor,
[WolfeSearch](https://statmodels7.github.io/optimizers7/reference/WolfeSearch-class.md)
and
[NonmonotoneSearch](https://statmodels7.github.io/optimizers7/reference/NonmonotoneSearch-class.md)
for the other two.
