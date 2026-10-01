# S7 Class for Chebyshev's Method

A
[Newton](https://statmodels7.github.io/optimizers7/reference/Newton-class.md)
optimizer whose direction carries Chebyshev's third-order correction.
Built by
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md).

## Usage

``` r
Chebyshev(
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
  typical = NULL,
  ratio = integer(0)
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

  The smallest eigenvalue the repaired Hessian may have.

- max_length:

  The largest component that a step may have.

- typical:

  `NULL`, or the typical size of each parameter, in whose units the step
  lengths are read.

- ratio:

  The largest \\\lVert d - d_N\rVert / \lVert d_N\rVert\\ at which the
  corrected direction \\d\\ is used.

## Value

An S7 object of class `Chebyshev` inheriting from
[Newton](https://statmodels7.github.io/optimizers7/reference/Newton-class.md).

## Details

Beyond the properties of
[Newton](https://statmodels7.github.io/optimizers7/reference/Newton-class.md),
a `Chebyshev` carries `ratio`, the largest relative length the
correction may have before the Newton step is kept instead.

## See also

[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
for the constructor.
