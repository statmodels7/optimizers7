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
  ncores = NULL,
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

- ncores:

  How many processes to spread the starts over. Defaults to `NULL`,
  meaning as many as are worth using: see Details. Pass `1` to stay in
  this session.

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

### Parallelism, which is arranged for you

The starts are independent, so they are run in parallel, and there is
nothing to set up: `ncores` is a number, and everything behind it —
starting the workers, loading the package on them, giving each an
independent random stream, and shutting them down again — happens inside
and is undone before the function returns. The default is
`min(n, max(1, parallel::detectCores() - 2))`: as many processes as
there are starts to run, but never so many that the machine has nothing
left for anything else. Two cores are held back rather than one because
the process doing the asking is one of them.

The mechanism differs by platform and the difference is not visible from
outside. On Linux and macOS the workers are forks, which start in
microseconds and inherit this session entire — including a package
loaded with pkgload, which is why development works there without
installing anything. Windows has no `fork`, so a socket cluster is
started instead and optimizers7 is loaded on each worker; if that fails,
because the package is not installed anywhere the workers can see it,
the run says so and continues sequentially rather than failing.

Threading inside C++ is not an option, and the reason is worth stating
because it is not about the objective. **The stopping rule is an R
object**, by design, and it is consulted at every iteration — that is
the feature this package was built around. Every run therefore returns
to R, and R is single-threaded, so calling it from several threads is
undefined behaviour that crashes rather than errors. Separate processes
have no such problem.

[`set.seed`](https://rdrr.io/r/base/Random.html) reproduces a parallel
run exactly, and reproduces it across platforms and across values of
`ncores`: the starting points are drawn here, before anything is
dispatched, and each worker is handed a stream derived from this
session's, so the same seed gives the same answer whether the work was
split eight ways or not at all.

### A start that fails is not a run that fails

A random start can easily land where the objective is undefined. Such a
start is recorded as failed and the others carry on; only if *every*
start fails is that an error. Otherwise one bad draw would throw away
nineteen good answers.

## See also

[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)

## Examples

``` r
multistart(bfgs())
#> <optimizer> multistart (BFGS)
#>   stop when : gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
#>   budgets   : maxit 10, evaluations 20000
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
