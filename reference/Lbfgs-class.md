# S7 Class for Limited-Memory BFGS

This class is an optimizer that holds the number of secant pairs to
keep, the curvature threshold that a pair must meet to be stored, and
the step length and the line search.
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
builds it. The pairs themselves exist only inside the run.

## Usage

``` r
Lbfgs(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  line_search = NULL,
  memory = integer(0),
  curv_tol = integer(0)
)
```

## Arguments

- step, line_search:

  The initial step length and the
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object, as in
  [`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md).

- memory:

  How many secant pairs to keep.

- curv_tol:

  The curvature threshold below which a pair is not stored.

## Value

An S7 object of class `Lbfgs` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the four properties above beside the seven shared ones.

## Details

Beyond the seven properties of every optimizer, an `Lbfgs` carries four:
`step` and `line_search`, and `memory` and `curv_tol`, which are its
own. It has no counterpart of `max_skip` of
[Bfgs](https://statmodels7.github.io/optimizers7/reference/Bfgs-class.md):
a pair that fails the curvature test is not stored (the trace records
`lbfgs update skipped`), and there is no accumulated matrix to protect.
When the direction is not a descent direction, the stored pairs are
cleared (`lbfgs reset`).

## See also

[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
for the constructor,
[Bfgs](https://statmodels7.github.io/optimizers7/reference/Bfgs-class.md)
for the full-matrix version.
