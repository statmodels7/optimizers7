# S7 Class for Nelder-Mead

An optimizer holding the size and shape of the initial simplex, which
set of reflection coefficients to use, and the safeguard against a
simplex that has collapsed. Built by
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md).
The simplex itself moves inside the run; the `simplex` property is a
starting one, or `NULL`.

## Usage

``` r
NelderMead(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  adaptive = logical(0),
  max_restarts = integer(0),
  degenerate_tol = integer(0),
  simplex = NULL
)
```

## Arguments

- step:

  Relative size of the initial simplex.

- adaptive:

  Logical; whether the dimension-dependent coefficients of Gao and Han
  are used.

- max_restarts:

  How many times a degenerate simplex may be rebuilt.

- degenerate_tol:

  The conditioning below which it is rebuilt.

- simplex:

  An optional starting simplex, a matrix with one vertex per row, or
  `NULL`.

## Value

An S7 object of class `NelderMead` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the five properties above beside the seven shared ones.

## Details

Beyond the seven properties every optimizer has, a `NelderMead` carries
five of its own: `step` and `simplex` describe where the run begins,
`adaptive` which coefficients it uses, and `max_restarts` with
`degenerate_tol` the safeguard.

## See also

[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
for the constructor,
[Compass](https://statmodels7.github.io/optimizers7/reference/Compass-class.md)
for the other derivative-free method.
