# Assemble the Result of a Chain

Combines the per-stage results into one: the last stage's point, value,
gradient, iteration count and verdict, the evaluation counts added over
every stage, the traces stacked with a `stage` column, and the messages
joined with their stage numbers.

## Usage

``` r
build_chain_result(results, optimizer, elapsed)
```

## Arguments

- results:

  The per-stage
  [`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
  objects, in order.

- optimizer:

  The `Chain` that produced them.

- elapsed:

  Total seconds for the whole chain.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

## Details

The counts are added with `Reduce("+")` rather than summed coordinate by
coordinate, so that the storage type survives: a chain of one stage must
report that stage's counts identically, and a
[`sum()`](https://rdrr.io/r/base/sum.html) over a
[`vapply()`](https://rdrr.io/r/base/lapply.html) would hand back doubles
where the stage had integers.

Traces are stacked only when every stage reports the same columns. A
derivative-free stage has a `stationarity` column where a descent has
`gnorm`, and when they disagree the last stage's trace is kept alone.
