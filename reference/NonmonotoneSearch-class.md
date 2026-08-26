# S7 Class for the Nonmonotone Line Search

Armijo backtracking whose reference value is the worst of the last
`memory + 1` objective values instead of the current one, so a step may
make things worse now to be better placed later. Built by
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md).

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

  The factor the step is multiplied by on each backtrack.

- memory:

  How many earlier values the reference looks back over.

- max_step:

  The most backtracks allowed. A **count**.

- resolution:

  What the objective can tell apart, a number or a function of no
  arguments.

## Value

An S7 object of class `NonmonotoneSearch` inheriting from
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md),
with the five properties above beside `label`.

## Details

Beyond the `label` every line search carries, a `NonmonotoneSearch`
holds `c1`, `shrink`, `memory`, `max_step` and `resolution`. `memory` is
the one property no other line search has, and at `memory = 0` the
object behaves exactly as an
[ArmijoSearch](https://statmodels7.github.io/optimizers7/reference/ArmijoSearch-class.md)
with the same constants.

## See also

[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
for the constructor,
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) for
the method that needs it.
