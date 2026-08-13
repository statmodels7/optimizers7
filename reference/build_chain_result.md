# Assemble the Result of a Chain

The last stage's point and verdict, with the work of every stage summed
and the traces stacked.

## Usage

``` r
build_chain_result(results, optimizer, elapsed)
```

## Arguments

- results:

  The per-stage results, in order.

- optimizer:

  The `Chain`.

- elapsed:

  Total seconds.

## Value

An
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).
