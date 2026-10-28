# S7 Class for Conjugate Gradients

A `Cg` object holds which of the four \\\beta\\ formulas combines the
gradient with the previous direction, how often the method is restarted
at steepest descent, and the step length and the line search. It is
built by
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
  line_search = NULL,
  stop_decrement = integer(0)
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

- stop_decrement:

  The largest decrease at which a stopping rule on a change may end the
  run as converged; see
  [`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md).

## Value

An S7 object of class `Cg` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the five properties above beside the seven shared ones.

## Details

Beyond the seven properties of every optimizer, a `Cg` carries five:
`step`, `line_search` and `stop_decrement`, shared with the other
line-search methods, and `beta` and `restart_every`, which are its own.
The previous direction is not a property; it exists only inside the run.

## See also

[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md) for
the constructor,
[GradientDescent](https://statmodels7.github.io/optimizers7/reference/GradientDescent-class.md)
for steepest descent.
