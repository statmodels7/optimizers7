# How Many Processes to Use

Turns `ncores = NULL` into a number: as many processes as there are
starts, but not more than the number of cores minus two.

## Usage

``` r
resolve_ncores(ncores, n)
```

## Arguments

- ncores:

  The value given by the caller, possibly `NULL`.

- n:

  The number of starts.

## Value

A single integer, at least one.

## Details

The rule is `min(n, max(1, detectCores() - 2))`. Two cores are left
free, one for the calling session and one for the rest of the machine,
which otherwise stops responding. Starting more processes than there are
starts costs the start-up time of the extra processes, which on Windows
is seconds, whereas a fork takes microseconds.

It is capped at two under `R CMD check`, which sets
`_R_CHECK_LIMIT_CORES_` and fails a package that ignores it, and it
falls back to one process where `detectCores()` cannot tell.
