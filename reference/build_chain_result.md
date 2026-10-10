# Assemble the Result of a Chain

Combines the results of the stages into one: the point, value, gradient,
iteration count and verdict of the last stage, the evaluation counts
added over every stage, the traces stacked with a `stage` column, and
the messages joined with their stage numbers.

## Usage

``` r
build_chain_result(results, optimizer, elapsed)
```

## Arguments

- results:

  The
  [`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
  objects of the stages, in order.

- optimizer:

  The `Chain` that produced them.

- elapsed:

  Total seconds for the whole chain.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

## Details

The counts are added with `Reduce("+")` and not summed coordinate by
coordinate, so that their storage type is kept: a chain of one stage
reports the counts of that stage identically, whereas a
[`sum()`](https://rdrr.io/r/base/sum.html) over a
[`vapply()`](https://rdrr.io/r/base/lapply.html) would return doubles
where the stage had integers.

Traces are stacked only when every stage reports the same columns. A
derivative-free stage has a `stationarity` column where a descent has
`gnorm`, and when they differ only the trace of the last stage is kept.
