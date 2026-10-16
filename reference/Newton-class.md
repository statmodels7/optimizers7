# S7 Class for Newton's Method

This class is an optimizer that holds the initial step length, the line
search, the two settings that determine how an indefinite Hessian is
repaired, and the two settings that bound and scale a step.
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
builds it. Among the shipped methods, only this one and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md),
which inherits from it, read the `he` argument of
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Usage

``` r
Newton(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  line_search = NULL,
  hessian_mod = character(0),
  floor = integer(0),
  max_length = integer(0),
  typical = NULL
)
```

## Arguments

- step, line_search:

  The initial step length and the
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object, as in
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

- hessian_mod:

  How an indefinite Hessian is repaired, `"eigen"` or `"ridge"`.

- floor:

  The lower bound of the repair; see
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

- max_length:

  The largest component that a step may have.

- typical:

  `NULL`, or the typical size of each parameter, in the units of which
  the step lengths are measured.

## Value

An S7 object of class `Newton` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the six properties above beside the seven shared ones.

## Details

Beyond the seven properties of every optimizer, a `Newton` carries six:
`step` and `line_search`, shared with the other line-search methods, and
`hessian_mod`, `floor`, `max_length` and `typical`, which are its own.

## See also

[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
for the constructor,
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
for the method that needs no Hessian.
