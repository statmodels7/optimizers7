# Assemble the Result of a Multi-Start Run

Picks the run with the lowest value, keeps everything that run carried,
and adds what the others found: the per-start table, the counts, and the
number of distinct optima.

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

The count of distinct optima is the reason to run this at all, so it is
computed here: the values reached are sorted and cut wherever
consecutive ones differ by more than `distinct_tol`. It is a statement
about the objective, not about the optimizer, and it is the one piece of
evidence a single run can never supply.
