# S7 Class for the Barzilai-Borwein Method

A `Bb` object holds which of the two Rayleigh quotients gives the step
length, the bounds to which that quotient is clamped, the curvature
threshold that a secant pair must meet, and the step multiplier and the
line search. It is built by
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md). Its
line search defaults to
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md);
it is the only shipped method with that default.

## Usage

``` r
Bb(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  variant = character(0),
  alpha0 = integer(0),
  alpha_min = integer(0),
  alpha_max = integer(0),
  curv_tol = integer(0),
  step = integer(0),
  line_search = NULL
)
```

## Arguments

- variant:

  Which step-length formula, one of `"alternate"`, `"bb1"`, `"bb2"`.

- alpha0:

  The step length used before there is a secant pair.

- alpha_min, alpha_max:

  The bounds to which the step length is clamped.

- curv_tol:

  The relative threshold below which a secant pair is rejected.

- step, line_search:

  The multiplier offered to the line search and the
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object, as in
  [`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md).

## Value

An S7 object of class `Bb` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the seven properties above beside the seven shared ones.

## Details

Beyond the seven properties of every optimizer, a `Bb` carries seven of
its own: `variant`, `alpha0`, `alpha_min`, `alpha_max` and `curv_tol`
for the step-length estimate, and `step` and `line_search` for the
search along the direction.

## See also

[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) for
the constructor,
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
for the acceptance test that is its default.
