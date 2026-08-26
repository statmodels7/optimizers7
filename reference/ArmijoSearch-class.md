# S7 Class for Armijo Backtracking

A line search that starts from the full step and halves it until the
objective falls by enough. Built by
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md).
It evaluates the objective at trial points and never the gradient, which
makes it the cheap choice per iteration.

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

  The factor the step is multiplied by on each backtrack.

- max_step:

  The most backtracks allowed. A **count**, not a length.

- resolution:

  What the objective can tell apart, a number or a function of no
  arguments.

## Value

An S7 object of class `ArmijoSearch` inheriting from
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md),
with the four properties above beside `label`.

## Details

Beyond the `label` every line search carries, an `ArmijoSearch` holds
`c1`, `shrink`, `max_step` and `resolution`. It has no `c2`, having no
curvature condition;
[`line_search_spec()`](https://statmodels7.github.io/optimizers7/reference/line_search_spec.md)
fills that field with a value the compiled side ignores so that all
three searches describe themselves in the same shape.

## See also

[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
for the constructor,
[WolfeSearch](https://statmodels7.github.io/optimizers7/reference/WolfeSearch-class.md)
and
[NonmonotoneSearch](https://statmodels7.github.io/optimizers7/reference/NonmonotoneSearch-class.md)
for the other two.
