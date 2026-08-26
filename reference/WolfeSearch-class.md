# S7 Class for the Strong Wolfe Line Search

A line search that brackets a step satisfying both Wolfe conditions and
then bisects inside the bracket. Built by
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md).
It evaluates the **gradient** at trial points as well as the objective.
That is what the curvature condition costs, and what a quasi-Newton
method needs.

## Usage

``` r
WolfeSearch(
  label = character(0),
  c1 = integer(0),
  c2 = integer(0),
  max_step = integer(0),
  resolution = NULL
)
```

## Arguments

- c1:

  The sufficient-decrease constant.

- c2:

  The curvature constant.

- max_step:

  The most trial steps allowed in each phase. A **count**.

- resolution:

  What the objective can tell apart, a number or a function of no
  arguments.

## Value

An S7 object of class `WolfeSearch` inheriting from
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md),
with the four properties above beside `label`.

## Details

Beyond the `label` every line search carries, a `WolfeSearch` holds
`c1`, `c2`, `max_step` and `resolution`. It has no `shrink`, taking no
backtracking steps of a fixed ratio;
[`line_search_spec()`](https://statmodels7.github.io/optimizers7/reference/line_search_spec.md)
fills that field with a value the compiled side ignores.

## See also

[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
for the constructor,
[ArmijoSearch](https://statmodels7.github.io/optimizers7/reference/ArmijoSearch-class.md)
and
[NonmonotoneSearch](https://statmodels7.github.io/optimizers7/reference/NonmonotoneSearch-class.md)
for the other two.
