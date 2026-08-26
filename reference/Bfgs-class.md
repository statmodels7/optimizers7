# S7 Class for BFGS

An optimizer holding the initial step length, the line search, and the
two settings that govern what happens when a secant pair carries no
usable curvature. Built by
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md).
The inverse-Hessian approximation itself is not a property: it is built
inside the run and discarded with it.

## Usage

``` r
Bfgs(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  line_search = NULL,
  curv_tol = integer(0),
  max_skip = integer(0)
)
```

## Arguments

- step, line_search:

  The initial step length and the
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object, as in
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md).

- curv_tol:

  The curvature threshold below which the update is skipped.

- max_skip:

  Consecutive skips before the approximation is reset to the identity.

## Value

An S7 object of class `Bfgs` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the four properties above beside the seven shared ones.

## Details

Beyond the seven properties every optimizer has, a `Bfgs` carries four:
`step` and `line_search`, shared with the other line-search methods, and
`curv_tol` and `max_skip`, which are its own.

## See also

[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
for the constructor,
[Lbfgs](https://statmodels7.github.io/optimizers7/reference/Lbfgs-class.md)
for the limited-memory version.
