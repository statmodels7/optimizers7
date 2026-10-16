# S7 Class for Pattern Search

A `Compass` object holds the poll size and how it changes, the set of
poll directions that is used, and whether the poll stops at the first
improvement. It is built by
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md).
With `directions = "mads"` a run draws from R's generator and records
its seed; with `"coordinate"` it is deterministic and records none.

## Usage

``` r
Compass(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0),
  step = integer(0),
  directions = character(0),
  opportunistic = logical(0),
  expand = integer(0),
  shrink = integer(0)
)
```

## Arguments

- step:

  Initial poll size, relative to the starting value.

- directions:

  Either `"mads"` or `"coordinate"`.

- opportunistic:

  Logical; whether the poll stops at the first improvement that it
  finds.

- expand, shrink:

  Factors applied to the poll size after a success and after a failure.
  After a success the poll size is capped at \\10^4\\ times its initial
  value.

## Value

An S7 object of class `Compass` inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
with the five properties above beside the seven shared ones.

## Details

Beyond the seven properties of every optimizer, a `Compass` carries five
of its own: `step`, `expand` and `shrink` govern the poll radius,
`directions` selects the set that is polled, and `opportunistic`
controls when the poll stops.

## See also

[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
for the constructor,
[NelderMead](https://statmodels7.github.io/optimizers7/reference/NelderMead-class.md)
for the other derivative-free method.
