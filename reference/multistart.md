# Run an Optimizer From Many Starting Points

Wraps any optimizer and runs it from several starting points, returning
the best result together with the number of distinct optimal values
reached. Two runs count as reaching the same optimum when their
objective values differ by at most `distinct_tol`. A count above one
shows that the surface has more than one local minimum; a count of one
does not exclude other minima with the same value or minima that all the
starts missed.

## Usage

``` r
multistart(
  optimizer,
  n = 10,
  starts = NULL,
  spread = 1,
  ncores = NULL,
  distinct_tol = 1e-06,
  verbose = FALSE,
  refresh = 1,
  keep_trace = TRUE
)
```

## Arguments

- optimizer:

  The
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  to run, any of them, including another `multistart()`. Anything else
  raises an error naming
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  as an example.

- n:

  How many starts, the caller's own `par` among them. Defaults to `10`.
  It becomes the `maxit` of the wrapper, so `n = 0` is rejected with an
  error that names `maxit`.

- starts:

  An optional matrix of starting points, one per row, used verbatim; `n`
  and `spread` are then ignored and `n` is taken from the number of
  rows. A matrix whose number of columns differs from the length of
  `par` is rejected when the run starts, with an error whose message
  gives the required number of columns.

- spread:

  How widely the random starts are scattered, in units of the
  unconstrained scale. Defaults to `1`.

- ncores:

  The number of processes over which to spread the starts. Defaults to
  `NULL`, meaning `min(n, max(1, parallel::detectCores() - 2))`. Pass
  `1` for a sequential run. For a deterministic inner optimizer the
  answer does not depend on the value; for a stochastic one, such as
  [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md),
  it does.

- distinct_tol:

  Objective values differing by less than this count as the same
  optimum. Defaults to `1e-6`.

- verbose:

  Report each start as it finishes? Defaults to `FALSE`. Only a
  sequential run reports; a run over several processes prints nothing.

- refresh:

  Report every this many starts. Defaults to `1`.

- keep_trace:

  Keep the per-start summary? Defaults to `TRUE`, unlike every other
  optimizer, because here the trace is one row per start, not one per
  iteration.

## Value

An S7 object of class
[MultiStart](https://statmodels7.github.io/optimizers7/reference/MultiStart-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## What the result reports

The result has the `par`, `value`, `gradient`, `converged` and
`criterion_met` of the best run. Its `counts` are the sums over all
successful starts, its `iterations` is the number of starts, and its
`message` and `trace` are those of the wrapper. The summary of the
starts is in `trace`, one row for each start, with the columns `start`,
`value`, `converged` and `iterations`. The `message` counts the starts
that succeeded and those that converged, the distinct optimal values
found and how often the best was reached, as in

    12 starts, 12 succeeded, 12 converged, 2 distinct optima;
    the best was found 6 times.

The count is computed by sorting the values reached and cutting wherever
consecutive ones differ by more than `distinct_tol`.

## Starting points

Unless `starts` is given, the first starting point is `par`, and the
remaining \\n - 1\\ form a Latin hypercube: the range of each coordinate
is divided into equal strata and each stratum is used exactly once,
which spreads the starts more evenly than independent draws. They are
generated on the unconstrained scale and mapped back through the bounds,
so every start is admissible by construction, however tight the box. The
rows of a matrix given in `starts` are not checked against the bounds in
advance: a row outside the box fails in its own run and is recorded as a
failed start.

## Parallel execution

The starts are independent and are run over `ncores` processes. The
function creates the workers, loads the package, assigns the random
streams and shuts the workers down internally: on Unix-alikes the
workers are forks, on Windows they form a socket cluster. If the workers
cannot load optimizers7, the run gives a warning and proceeds
sequentially. If optimizers7 was loaded from source with pkgload, the
starts are run sequentially on every platform, without a warning,
because a worker would load an installed copy of the package instead of
the copy in this session.

Processes are used because the stopping rule is an R object evaluated at
every iteration, and R cannot be called from several threads.

The starting points are drawn in the calling session before dispatch,
and each worker receives a random stream derived from the seed of the
session, so [`set.seed()`](https://rdrr.io/r/base/Random.html)
reproduces the run. When the inner optimizer is deterministic, such as
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
or
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
the result is the same for every value of `ncores`. When the inner
optimizer draws random numbers, such as
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md), the
result depends on the assignment of the streams to the workers, and
therefore on `ncores`.

## Failed starts

A start where the objective is undefined is recorded as failed and the
others proceed; the `message` then reports how many succeeded and quotes
the first failure. An error is raised only when every start fails, and
it quotes the first error.

## References

McKay, M. D., Beckman, R. J. and Conover, W. J. (1979). A comparison of
three methods for selecting values of input variables in the analysis of
output from a computer code. *Technometrics* **21**, 239–245.

## See also

[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
for running optimizers in sequence instead of in parallel,
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) for
a global search that needs no restarts,
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md)
for another way of generating starting points.

## Examples

``` r
multistart(bfgs())
#> <optimizer> multistart (BFGS)
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 10, evaluations Inf
#>   settings  : optimizer = BFGS, n = 10, starts = <NULL>, spread = 1, ncores = <NULL>, distinct_tol = 1e-06

# The surface has two minima, and one of them is lower. A single run from the
# origin finds the lower one; the count shows that there is another.
f <- function(p) (p[1]^2 - 1)^2 + p[2]^2 + 0.3 * p[1]
set.seed(1)
r <- minimize(multistart(bfgs(), n = 12), f, c(0, 0))
r@message
#> [1] "12 starts, 12 succeeded, 12 converged, 2 distinct optima; the best was found 6 times. Best run: gradient obtained by finite differences"
table(round(r@trace$value, 6))
#> 
#> -0.305428  0.294146 
#>         6         6 
r@par
#> [1] -1.035579  0.000000

# With a deterministic inner optimizer the answer does not depend on the
# number of processes.
set.seed(11); one  <- minimize(multistart(bfgs(), n = 8, ncores = 1), f, c(0, 0))
set.seed(11); many <- minimize(multistart(bfgs(), n = 8, ncores = 2), f, c(0, 0))
identical(one@value, many@value)
#> [1] TRUE

# A start where the objective is undefined is recorded, not fatal.
g <- function(p) if (p[1] > 0.5) stop("undefined here") else sum(p^2)
set.seed(3)
minimize(multistart(bfgs(), n = 8, ncores = 1), g, c(0, 0))@message
#> [1] "8 starts, 5 succeeded, 5 converged, 1 distinct optimum; the best was found 5 times. First failure: undefined here. Best run: gradient obtained by finite differences"

# If every start fails, the run raises an error that quotes the first
# failure.
try(minimize(multistart(bfgs(), n = 4, ncores = 1),
             function(p) stop("never works"), c(0, 0)))
#> Error : Every start failed. The first error was: never works

# Starting points may be given outright, and then n is the row count.
S <- rbind(c(-2, 0), c(2, 0), c(0, 3))
minimize(multistart(bfgs(), starts = S, ncores = 1), f, c(0, 0))@trace
#>   start      value converged iterations
#> 1     1 -0.3054285      TRUE          5
#> 2     2  0.2941465      TRUE          6
#> 3     3 -0.3054285      TRUE          9
```
