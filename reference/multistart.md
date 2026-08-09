# Run an Optimizer From Many Starting Points

Wraps any optimizer and runs it from several starts, returning the best
result together with the number of distinct answers found.

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

  The optimizer to run. Any of them, including another `multistart()`.

- n:

  How many starts, the user's own `par` among them. Defaults to `10`.

- starts:

  An optional matrix of starting points, one per row, used verbatim; `n`
  and `spread` are then ignored.

- spread:

  How widely the random starts are scattered, in units of the
  unconstrained scale. Defaults to `1`.

- ncores:

  How many processes to spread the starts over. Defaults to `NULL`,
  meaning `min(n, parallel::detectCores() - 2)`. Pass `1` for a
  sequential run.

- distinct_tol:

  Objective values differing by less than this are counted as the same
  optimum. Defaults to `1e-6`.

- verbose:

  Report each start as it finishes? Defaults to `FALSE`.

- refresh:

  Report every this many starts. Defaults to `1`.

- keep_trace:

  Keep the per-start summary? Defaults to `TRUE` — it is one row per
  start, not one per iteration.

## Value

An S7 object of class `MultiStart`, inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Details

Every method in this package finds a local minimum. Running from several
starting points is the general remedy, and beyond the best point found
it reports the **count of distinct optima**, which is evidence about
whether the objective has a single minimum at all.

The result is the best run, with everything it carried. The per-start
summary is in `trace`: one row each, with the value reached, whether
that start converged, and how many iterations it took. The message
counts the starts that succeeded, the ones that converged, and the
distinct optima found.

### Starting points

The first starting point is always `par`. The remaining `n - 1` form a
Latin hypercube: each coordinate's range is divided into equal strata
and each stratum is used exactly once, which spreads the starts more
evenly than independent draws. They are generated on the *unconstrained*
scale and mapped back through the bounds, so every start is admissible
by construction.

### Parallel execution

The starts are independent and are run in parallel over `ncores`
processes. The default is `min(n, max(1, parallel::detectCores() - 2))`.
Worker creation, package loading, random-stream assignment and shutdown
are handled internally: on Unix-alikes the workers are forks, on Windows
a socket cluster, and if the workers cannot load the package the run
warns and proceeds sequentially. Processes are used rather than threads
because the stopping rule is an R object consulted at every iteration,
and R cannot be called from multiple threads.

The starting points are drawn in the calling session before dispatch,
and each worker receives a random stream derived from the session's
seed, so [`set.seed`](https://rdrr.io/r/base/Random.html) reproduces the
run identically for any value of `ncores` and on any platform.

### Failed starts

A start where the objective is undefined is recorded as failed and the
remaining starts proceed; an error is raised only when every start
fails.

## References

McKay, M. D., Beckman, R. J. and Conover, W. J. (1979). A comparison of
three methods for selecting values of input variables in the analysis of
output from a computer code. *Technometrics* **21**, 239–245.

## See also

[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)

## Examples

``` r
multistart(bfgs())
#> <optimizer> multistart (BFGS)
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-12 (relative)
#>   budgets   : maxit 10, evaluations Inf
#>   settings  : optimizer = BFGS, n = 10, starts = <NULL>, spread = 1, ncores = <NULL>, distinct_tol = 1e-06

# a surface with two minima, one of them better
f <- function(p) (p[1]^2 - 1)^2 + p[2]^2 + 0.3 * p[1]
set.seed(1)
r <- minimize(multistart(bfgs(), n = 12), f, c(0, 0))
r@message
#> [1] "12 starts, 12 succeeded, 12 converged, 2 distinct optima; the best was found 6 times. Best run: gradient obtained by finite differences"
table(round(r@trace$value, 6))
#> 
#> -0.305428  0.294146 
#>         6         6 
```
