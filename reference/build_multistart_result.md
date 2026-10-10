# Assemble the Result of a Multi-Start Run

Picks the run with the lowest value and takes its `par`, `value`,
`gradient`, `converged` and `criterion_met`. The evaluation counts are
summed over the successful starts, `iterations` is the number of starts,
and the table of the starts, the counts of starts and the number of
distinct optimal values are added in `trace` and `message`.

## Usage

``` r
build_multistart_result(res, S, optimizer, seed, elapsed)
```

## Arguments

- res:

  The list of results, entries that failed being character messages.

- S:

  The matrix of starting points.

- optimizer:

  The `MultiStart` object.

- seed:

  The generator state the run began with.

- elapsed:

  Seconds.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

## Details

The number of distinct optimal values is computed by sorting the values
reached and cutting wherever consecutive values differ by more than
`distinct_tol`. It describes the objective and not the optimizer, and a
single run cannot provide it.
