# S7 Class for BFGS

This class is an optimizer that holds the initial step length, the line
search, and the two settings that apply when a secant pair has no usable
curvature.
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
builds it. The inverse-Hessian approximation is not a property: it is
built inside the run and discarded at its end.

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
  max_skip = integer(0),
  stop_decrement = integer(0)
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

- stop_decrement:

  The largest decrease at which a stopping rule on a change may end the
  run as converged; see
  [`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md).

## Value

An S7 object of class `Bfgs` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the five properties above beside the seven shared ones.

## Details

Beyond the seven properties of every optimizer, a `Bfgs` carries five:
`step`, `line_search` and `stop_decrement`, shared with the other
line-search methods, and `curv_tol` and `max_skip`, which are its own.

## See also

[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
for the constructor,
[Lbfgs](https://statmodels7.github.io/optimizers7/reference/Lbfgs-class.md)
for the limited-memory version.
