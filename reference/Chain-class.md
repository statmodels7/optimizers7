# S7 Class for a Sequence of Optimizers

An optimizer holding a list of optimizers to be run in order, each from
the point the previous one reached. Built by
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md).
Because it inherits from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
it is accepted anywhere a single method is, including inside
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md).

## Usage

``` r
Chain(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  stages = list()
)
```

## Arguments

- stages:

  A list of
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  objects, in the order they run.

## Value

An S7 object of class `Chain` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the property `stages` beside the seven shared ones.

## Details

Beyond the seven properties every optimizer has, a `Chain` carries
`stages`, the list of optimizers in running order. The seven it inherits
describe the **last** stage: its `criterion`, `maxit` and `max_eval` are
copied from there, because the last stage is the one whose rule ends the
run and whose result is reported. `refresh` is fixed at 1, the chain's
own progress being one line per stage.

## See also

[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
for the constructor,
[`minimize.Chain()`](https://statmodels7.github.io/optimizers7/reference/minimize.Chain.md)
for the run.
