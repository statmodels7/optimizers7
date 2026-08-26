# S7 Class for Conjugate Gradients

An optimizer holding which of the four \\\beta\\ formulas bends the
direction, how often the method is restarted at steepest descent, and
the usual step length and line search. Built by
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md).

## Usage

``` r
Cg(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  beta = character(0),
  restart_every = integer(0),
  step = integer(0),
  line_search = NULL
)
```

## Arguments

- beta:

  Which update formula, one of `"pr"`, `"fr"`, `"hs"`, `"dy"`.

- restart_every:

  How often the method is restarted at steepest descent; `0` means
  never.

- step:

  The initial step length offered to the line search.

- line_search:

  A
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object.

## Value

An S7 object of class `Cg` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the four properties above beside the seven shared ones.

## Details

Beyond the seven properties every optimizer has, a `Cg` carries four:
`step` and `line_search`, shared with the other line-search methods, and
`beta` and `restart_every`, which are its own. The previous direction is
not a property; it lives in the run.

## See also

[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md) for
the constructor,
[GradientDescent](https://statmodels7.github.io/optimizers7/reference/GradientDescent-class.md)
for the method it improves on.
