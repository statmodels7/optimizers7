# S7 Class for Brent's Method

A `Brent` object is an optimizer for a function of one variable. It
holds the length of the first step of the bracketing phase and the
number of times that phase may expand, and it is built by
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md).

## Usage

``` r
Brent(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  max_expand = integer(0)
)
```

## Arguments

- step:

  The length of the first trial step.

- max_expand:

  The number of expansions allowed while bracketing.

## Value

An S7 object of class `Brent` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the two properties above beside the seven shared ones.

## Details

Beyond the seven properties of every optimizer, a `Brent` carries two of
its own: `step`, the length of the first trial step, and `max_expand`,
the number of golden expansions that the bracketing phase may take.

## See also

[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
for the constructor.
