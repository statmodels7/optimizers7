# S7 Class for Chebyshev's Method

This class is a
[Newton](https://statmodels7.github.io/optimizers7/reference/Newton-class.md)
optimizer whose direction carries Chebyshev's third-order correction.
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
builds it.

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
  modify = NULL,
  max_length = integer(0),
  typical_size = NULL,
  stop_decrement = integer(0),
  ratio = integer(0)
)
```

## Arguments

- step, line_search:

  The initial step length and the
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object, as in
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

- modify:

  The
  [`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md)
  object that repairs an indefinite Hessian; see
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

- max_length:

  The largest component that a step may have.

- typical_size:

  `NULL`, or the typical size of each parameter, in the units of which
  the step lengths are measured.

- stop_decrement:

  The largest Newton decrement at which a stopping rule on a change may
  end the run; see
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

- ratio:

  The largest \\\lVert d - d_N\rVert / \lVert d_N\rVert\\ at which the
  corrected direction \\d\\ is used.

## Value

An S7 object of class `Chebyshev` inheriting from
[Newton](https://statmodels7.github.io/optimizers7/reference/Newton-class.md).

## Details

Beyond the properties of
[Newton](https://statmodels7.github.io/optimizers7/reference/Newton-class.md),
a `Chebyshev` carries `ratio`, the largest relative length of the
correction at which the corrected direction is used; above it the Newton
step is used.

## See also

[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
for the constructor.
