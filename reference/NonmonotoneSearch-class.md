# S7 Class for the Nonmonotone Line Search

This class is Armijo backtracking whose reference value is the largest
of the last `memory + 1` objective values, including the current one, so
a step may increase the objective.
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
builds it.

## Usage

``` r
NonmonotoneSearch(
  label = character(0),
  c1 = integer(0),
  shrink = integer(0),
  memory = integer(0),
  max_step = integer(0),
  resolution = NULL
)
```

## Arguments

- c1:

  The sufficient-decrease constant.

- shrink:

  The factor by which the step is multiplied on each backtrack.

- memory:

  The number of earlier values, besides the current one, over which the
  reference is taken.

- max_step:

  The most trial points allowed, including the full step. It is a count.

- resolution:

  The smallest difference in the objective that the objective can
  resolve, a number or a function of no arguments.

## Value

An S7 object of class `NonmonotoneSearch` inheriting from
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md),
with the five properties above beside `label`.

## Details

Beyond the `label` of every line search, a `NonmonotoneSearch` holds
`c1`, `shrink`, `memory`, `max_step` and `resolution`. Only this class
has the property `memory`, and at `memory = 0` the object behaves
exactly as an
[ArmijoSearch](https://statmodels7.github.io/optimizers7/reference/ArmijoSearch-class.md)
with the same constants.

## See also

[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
for the constructor,
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) for
the method that uses it by default.
