# S7 Class for Multi-Start

An optimizer wrapping another and running it from several starting
points. Built by
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md).
Its seven shared properties describe the wrapper and not the inner run:
`criterion`, `max_eval` and the rest are copied from the optimizer
inside so that printing tells the truth, while `maxit` counts
**starts**.

## Usage

``` r
MultiStart(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  optimizer = NULL,
  n = integer(0),
  starts = NULL,
  spread = integer(0),
  ncores = NULL,
  distinct_tol = integer(0)
)
```

## Arguments

- optimizer:

  The inner
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
  run from each starting point.

- n:

  How many starts.

- starts:

  An optional matrix of starting points, one per row.

- spread:

  How widely the random starts are scattered.

- ncores:

  How many processes the starts are spread over.

- distinct_tol:

  Objective values closer than this count as one optimum.

## Value

An S7 object of class `MultiStart` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the six properties above beside the seven shared ones.

## Details

Beyond the seven shared properties a `MultiStart` carries six of its
own: `optimizer` is the one being wrapped, `n`, `starts` and `spread`
say where the runs begin, `ncores` how they are spread over processes,
and `distinct_tol` how close two answers must be to count as one.

Because the rule that is evaluated belongs to the inner optimizer,
[`with_criterion()`](https://statmodels7.github.io/optimizers7/reference/with_criterion.md)
has a method for this class that sets both; setting the outer property
alone changes the printing and nothing else.

## See also

[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
for the constructor,
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
for the other wrapper,
[`with_criterion()`](https://statmodels7.github.io/optimizers7/reference/with_criterion.md)
for why setting the outer rule is not enough.
