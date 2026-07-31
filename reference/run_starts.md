# Run the Starts, in Parallel or Not

Evaluates `one(i)` for `i` in `1:n`, over `ncores` processes, and cleans
up after itself.

## Usage

``` r
run_starts(one, n, ncores, verbose, refresh)
```

## Arguments

- one:

  A function of the start's index, returning a result or a message.

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

Three routes, chosen for the caller rather than by them.

One process is the sequential loop, and it is the only route that can
report progress as it goes, a worker having nowhere to print to that the
caller would see.

On a Unix-alike the workers are **forks**, through
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html). A
fork starts in microseconds and inherits this session entire, so there
is nothing to load and nothing to export — including a package loaded
with pkgload, which is why this works during development where a socket
cluster does not.

On Windows there is no `fork`, so a **socket cluster** is started here
and stopped on exit. Its workers are fresh sessions that know nothing,
so optimizers7 has to be loaded on them; when it cannot be, because it
is not installed anywhere they can see, this warns and runs sequentially
rather than failing. A slower answer beats an error about
`checkForRemoteErrors` for someone who only asked for several starting
points.

Both parallel routes seed their workers from this session's stream, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces the run
and reproduces it identically whatever `ncores` was. The generator is
set to `"L'Ecuyer-CMRG"` for the duration and put back afterwards, that
being the only kind R can split into independent streams.
