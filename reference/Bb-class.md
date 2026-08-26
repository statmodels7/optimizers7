# S7 Class for the Barzilai-Borwein Method

An optimizer holding which of the two Rayleigh quotients gives the step
length, the bounds that quotient is clamped to, the curvature threshold
a secant pair must meet, and the usual step multiplier and line search.
Built by
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md). Its
line search defaults to
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md),
the only shipped method for which it does.

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

  Bounds the step length is clamped to.

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

Beyond the seven properties every optimizer has, a `Bb` carries seven of
its own: `variant`, `alpha0`, `alpha_min`, `alpha_max` and `curv_tol`
for the step-length estimate, and `step` and `line_search` for what is
done with it.

## See also

[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) for
the constructor,
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
for the acceptance test it defaults to.
