# Run the Starts, in Parallel or Not

Evaluates `one(i)` for `i` in `1:n`, over `ncores` processes, and cleans
up after itself.

## Usage

``` r
run_starts(one, n, ncores, verbose, refresh)
```

## Arguments

- one:

  A function of the index of the start, returning a result or a message.

- n:

  How many starts.

- ncores:

  How many processes, already resolved.

- verbose:

  Report each start as it finishes?

- refresh:

  Report every this many.

## Value

A list of `n` results.

## Details

The function selects one of three routes, so that the caller does not
have to.

One process is the sequential loop, and it is the only route that
reports progress during the run, since a worker has no output that the
caller would see. When optimizers7 was loaded with pkgload, this route
is used whatever `ncores` is.

On a Unix-alike the workers are forks, through
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html). A
fork starts in microseconds and inherits the whole session, so the
function does not have to load the package or export objects.

On Windows there is no `fork`, so a socket cluster is started here and
stopped on exit. Its workers are new sessions, so optimizers7 has to be
loaded on them; when that is not possible, because the package is not
installed in a library that they can see, the function gives a warning
and runs the starts sequentially instead of failing.

Both parallel routes seed their workers from the stream of this session,
so [`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces the
run. For a deterministic inner optimizer the result is the same for
every `ncores`; for one that draws random numbers, such as
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md), it
depends on the assignment of the streams to the workers. The generator
is set to `"L'Ecuyer-CMRG"` for the duration and restored afterwards,
since it is the only kind that R can split into independent streams.
