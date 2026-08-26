# Run an Optimizer From Many Starting Points

Wraps any optimizer and runs it from several starting points, returning
the best result together with the number of **distinct optima** found.
That count is evidence about the objective that no single run can
supply: it says whether the surface has one minimum or several.

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
  It becomes the wrapper's `maxit`, so `n = 0` is refused with a message
  naming `maxit`.

- starts:

  An optional matrix of starting points, one per row, used verbatim; `n`
  and `spread` are then ignored and `n` is taken from the number of
  rows. A matrix whose column count does not match the parameter vector
  is refused when the run starts, naming both.

- spread:

  How widely the random starts are scattered, in units of the
  unconstrained scale. Defaults to `1`.

- ncores:

  How many processes to spread the starts over. Defaults to `NULL`,
  meaning `min(n, max(1, parallel::detectCores() - 2))`. Pass `1` for a
  sequential run; the answer does not depend on the value.

- distinct_tol:

  Objective values differing by less than this count as the same
  optimum. Defaults to `1e-6`.

- verbose:

  Report each start as it finishes? Defaults to `FALSE`.

- refresh:

  Report every this many starts. Defaults to `1`.

- keep_trace:

  Keep the per-start summary? Defaults to `TRUE`, unlike every other
  optimizer, because here the trace is one row per **start**, not one
  per iteration.

## Value

An S7 object of class
[MultiStart](https://statmodels7.github.io/optimizers7/reference/MultiStart-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## What the result reports

The result is the best run, with everything that run carried. The
per-start summary is in `trace`, one row each with the columns `start`,
`value`, `converged` and `iterations`. The `message` counts the starts
that succeeded, the ones that converged, the distinct optima found and
how often the best was reached, as in

    12 starts, 12 succeeded, 12 converged, 2 distinct optima;
    the best was found 6 times.

The count is computed by sorting the values reached and cutting wherever
consecutive ones differ by more than `distinct_tol`.

## Starting points

The first starting point is always `par`. The remaining \\n - 1\\ form a
Latin hypercube: each coordinate's range is divided into equal strata
and each stratum is used exactly once, which spreads the starts more
evenly than independent draws would. They are generated on the
**unconstrained** scale and mapped back through the bounds, so every
start is admissible by construction however tight the box.

## Parallel execution

The starts are independent and are run over `ncores` processes. Worker
creation, package loading, random-stream assignment and shutdown are
handled internally: on Unix-alikes the workers are forks, on Windows a
socket cluster. If the workers cannot load the package the run warns and
proceeds sequentially. That is what happens under `pkgload`, the workers
being separate sessions that do not inherit this one's loaded packages.

Processes are used because the stopping rule is an R object consulted at
every iteration, and R cannot be called from several threads.

The starting points are drawn in the calling session before dispatch,
and each worker receives a random stream derived from the session's
seed, so [`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces
the run identically at any `ncores`: measured on the example below,
`ncores` of 1, 2 and 4 give `-0.3054284837` and the same distinct count.

## Failed starts

A start where the objective is undefined is recorded as failed and the
rest proceed; the `message` then reports how many succeeded and quotes
the first failure. An error is raised only when **every** start fails,
and it quotes the first message.

## References

McKay, M. D., Beckman, R. J. and Conover, W. J. (1979). A comparison of
three methods for selecting values of input variables in the analysis of
output from a computer code. *Technometrics* **21**, 239–245.

## See also

[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
for running optimizers in sequence rather than in parallel,
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) for
a global search that needs no restarts,
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md)
for the other way of generating starting points.

## Examples

``` r
multistart(bfgs())
#> <optimizer> multistart (BFGS)
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 10, evaluations Inf
#>   settings  : optimizer = BFGS, n = 10, starts = <NULL>, spread = 1, ncores = <NULL>, distinct_tol = 1e-06

# A surface with two minima, one of them better. A single run from the
# origin happens to find the good one; the count is what says there is
# another.
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

# The answer does not depend on how many processes ran it.
set.seed(11); one  <- minimize(multistart(bfgs(), n = 8, ncores = 1), f, c(0, 0))
set.seed(11); many <- minimize(multistart(bfgs(), n = 8, ncores = 2), f, c(0, 0))
identical(one@value, many@value)
#> [1] TRUE

# A start where the objective is undefined is recorded, not fatal.
g <- function(p) if (p[1] > 0.5) stop("undefined here") else sum(p^2)
set.seed(3)
minimize(multistart(bfgs(), n = 8, ncores = 1), g, c(0, 0))@message
#> [1] "8 starts, 5 succeeded, 5 converged, 1 distinct optimum; the best was found 5 times. First failure: undefined here. Best run: gradient obtained by finite differences"

# Every start failing is fatal, and the message says why.
try(minimize(multistart(bfgs(), n = 4, ncores = 1),
             function(p) stop("never works"), c(0, 0)))
#> Error : Every start failed. The first said: never works

# Starting points may be given outright, and then n is the row count.
S <- rbind(c(-2, 0), c(2, 0), c(0, 3))
minimize(multistart(bfgs(), starts = S, ncores = 1), f, c(0, 0))@trace
#>   start      value converged iterations
#> 1     1 -0.3054285      TRUE          5
#> 2     2  0.2941465      TRUE          6
#> 3     3 -0.3054285      TRUE          9
```
