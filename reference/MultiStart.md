# Run an Optimiser From Many Starting Points

Wraps any optimiser and runs it from several starts, returning the best
result — and, as importantly, telling you how many different answers it
found.

## Usage

``` r
multistart(
  optimizer,
  n = 10,
  starts = NULL,
  spread = 1,
  cluster = NULL,
  distinct_tol = 1e-06,
  verbose = FALSE,
  refresh = 1,
  keep_trace = TRUE
)
```

## Arguments

- optimizer:

  The optimiser to run. Any of them, including another `multistart()`.

- n:

  How many starts, the user's own `par` among them. Defaults to `10`.

- starts:

  An optional matrix of starting points, one per row, used verbatim; `n`
  and `spread` are then ignored.

- spread:

  How widely the random starts are scattered, in units of the
  unconstrained scale. Defaults to `1`.

- cluster:

  An optional cluster from parallel. Defaults to `NULL`, meaning run the
  starts one after another; see Details.

- distinct_tol:

  Objective values differing by less than this are counted as the same
  optimum. Defaults to `1e-6`.

- verbose:

  Report each start as it finishes? Defaults to `FALSE`.

- refresh:

  Report every this many starts. Defaults to `1`.

- keep_trace:

  Keep the per-start summary? Defaults to `TRUE` — it is one row per
  start, not one per iteration, and it is the point of the method.

## Value

An S7 object of class `MultiStart`, inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Details

Every method in this package finds a *local* minimum, and none of them
can tell you whether it is the global one; that is not a shortcoming of
the implementations but of local optimisation. Running from several
starts is the only general answer, and the most valuable thing it
produces is not the best point but the **count of distinct optima**,
which is direct evidence about whether the question you asked has one
answer. A run reporting one optimum from twenty starts is worth far more
than a single run reporting convergence.

The result is the best run, with everything it carried. The per-start
summary is in `trace`: one row each, with the value reached, whether
that start converged, and how many iterations it took. The message
counts the starts that succeeded, the ones that converged, and the
distinct optima found.

### Where the starts come from

The first is always `par`: your guess is a hypothesis worth testing, and
silently discarding it would be rude.

The rest are a Latin hypercube — each coordinate's range is cut into `n`
equal strata and each is used exactly once, so the starts cannot all
cluster in one corner the way independent draws can. They are generated
on the *unconstrained* scale and mapped back, which is what makes bounds
automatic: a start for a variance is drawn as a log and comes back
positive, and no draw is ever rejected for being outside the box.

### Parallelism, and why it is at the level of processes

Pass a cluster made by
[`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)
and the starts are distributed across it. optimizers7 is loaded on the
workers for you; a cluster you made is a cluster you keep, so stopping
it is yours to do.

Threading the starts inside C++ is not possible, and the reason is worth
stating because it is not about the objective. **The stopping rule is an
R object**, by design, and it is consulted at every iteration — that is
the feature this package was built around, a criterion the user can
write. Every run therefore returns to R, and R is single-threaded, so
calling it from several threads is undefined behaviour that crashes
rather than errors.

Separate processes have no such problem, so that is the route taken.

For reproducibility across workers use
[`parallel::clusterSetRNGStream()`](https://rdrr.io/r/parallel/RngStream.html);
ordinary [`set.seed()`](https://rdrr.io/r/base/Random.html) governs only
the sequential path. The `seed` recorded in the result is the master's
state, and reproduces a *sequential* run only: the workers draw from
streams of their own, which the master never sees. Reproducing a cluster
run means recording the argument you gave `clusterSetRNGStream()`, and
there is no way for this function to do that for you.

### A start that fails is not a run that fails

A random start can easily land where the objective is undefined. Such a
start is recorded as failed and the others carry on; only if *every*
start fails is that an error. Otherwise one bad draw would throw away
nineteen good answers.

## See also

[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/Bfgs.md)

## Examples

``` r
multistart(bfgs())
#> <optimizer> multistart (BFGS)
#>   stop when : gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
#>   budgets   : maxit 10, evaluations 20000
#>   settings  : optimizer = BFGS, n = 10, starts = <NULL>, spread = 1, cluster = <NULL>, distinct_tol = 1e-06

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
