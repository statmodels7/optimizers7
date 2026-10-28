# S7 Class for Gradient Descent

A `GradientDescent` object holds the initial step length and the line
search only, since steepest descent does not keep a model of the
surface. It is built by
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
  line_search = NULL,
  stop_decrement = integer(0)
)
```

## Arguments

- step:

  The initial step length offered to the line search.

- line_search:

  A
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object.

- stop_decrement:

  The largest decrease at which a stopping rule on a change may end the
  run as converged; see
  [`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md).

## Value

An S7 object of class `GradientDescent` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with `step`, `line_search` and `stop_decrement` beside the seven shared
properties.

## Details

Beyond the seven properties of every optimizer, a `GradientDescent`
carries `step`, `line_search` and `stop_decrement`, shared with the
other line-search methods.

## See also

[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) for
the constructor,
[Cg](https://statmodels7.github.io/optimizers7/reference/Cg-class.md)
and
[Bb](https://statmodels7.github.io/optimizers7/reference/Bb-class.md)
for the two other first-order methods.
