# How Many Processes to Use

Turns `ncores = NULL` into a number: as many processes as there are
starts, but never more than the machine can spare.

## Usage

``` r
resolve_ncores(ncores, n)
```

## Arguments

- ncores:

  What the caller asked for, possibly `NULL`.

- n:

  The number of starts.

## Value

A single integer, at least one.

## Details

The rule is `min(n, max(1, detectCores() - 2))`. Two are held back
rather than one because the session doing the asking is itself one of
them, and a machine with nothing left over is a machine that stops
responding. Asking for more processes than there are starts wastes the
cost of starting them, which on Windows is seconds rather than
microseconds.

It is capped at two under `R CMD check`, which sets
`_R_CHECK_LIMIT_CORES_` and fails a package that ignores it, and it
falls back to one process where `detectCores()` cannot tell.
