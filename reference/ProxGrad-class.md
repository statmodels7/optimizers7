# S7 Class for the Proximal Gradient Method

An optimizer holding the two descriptions of the non-smooth part of an
objective, its proximal operator and its value, together with the three
settings of the accelerated iteration. Built by
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md).
It is the one shipped class whose
[`optimizer_bounded()`](https://statmodels7.github.io/optimizers7/reference/optimizer_bounded.md)
is `FALSE`: its constraint travels inside `prox`, so box bounds are
refused.

## Usage

``` r
ProxGrad(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  prox = function() NULL,
  g = function() NULL,
  accelerate = logical(0),
  step = integer(0),
  shrink = integer(0),
  restart = logical(0)
)
```

## Arguments

- prox:

  The proximal operator of the non-smooth part, `prox(v, step)`.

- g:

  The value of the non-smooth part, `g(par)`.

- accelerate:

  Logical; whether the momentum extrapolation is applied.

- step:

  The initial step length offered to the backtracking search.

- shrink:

  The factor a rejected step is multiplied by.

- restart:

  Logical; whether an increase in the objective resets the momentum.

## Value

An S7 object of class `ProxGrad` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the six properties above beside the seven shared ones.

## Details

Beyond the seven properties every optimizer has, a `ProxGrad` carries
six of its own. `prox` and `g` describe the same non-smooth term from
two sides and are both required. `accelerate`, `step`, `shrink` and
`restart` govern the iteration.

## See also

[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
for the constructor,
[`minimize.ProxGrad()`](https://statmodels7.github.io/optimizers7/reference/minimize.ProxGrad.md)
for the run.
