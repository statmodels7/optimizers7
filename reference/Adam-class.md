# S7 Class for Adam

An optimizer holding the learning rate, the two moment decay rates and
the three safeguards of the Adam iteration. Built by
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md).
It is the one shipped method whose default stopping rule is
[`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md),
so a run ends on its iteration budget and reports `converged = FALSE`.

## Usage

``` r
Adam(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  alpha = integer(0),
  beta1 = integer(0),
  beta2 = integer(0),
  eps = integer(0),
  decay = integer(0),
  amsgrad = logical(0)
)
```

## Arguments

- alpha:

  The learning rate, the size of a step when the gradient is steady.

- beta1, beta2:

  Decay rates for the first and second moment estimates.

- eps:

  Added to the square-rooted second moment before dividing.

- decay:

  Rate at which the learning rate is reduced.

- amsgrad:

  Logical; whether the second moment is held at its running maximum.

## Value

An S7 object of class `Adam` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the six properties above beside the seven shared ones.

## Details

Beyond the seven properties every optimizer has, an `Adam` carries six
of its own: `alpha`, `beta1` and `beta2` for the iteration as Kingma and
Ba published it, and `eps`, `decay` and `amsgrad` for the three repairs.

## See also

[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
for the constructor,
[`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md)
for its default rule.
