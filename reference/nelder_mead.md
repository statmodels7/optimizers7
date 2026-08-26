# The Nelder-Mead Simplex Method

The Nelder-Mead simplex method: \\p+1\\ points are maintained and the
worst is reflected, expanded, contracted or shrunk through the centroid
of the others. Only objective values are used; no derivative is
required.

## Usage

``` r
nelder_mead(
  criterion = crit_stationary(),
  step = 0.1,
  adaptive = TRUE,
  max_restarts = 3,
  degenerate_tol = 1e-06,
  simplex = NULL,
  maxit = 2000,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 50,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object. Defaults to `crit_stationary(1e-8)`, read on the diameter of
  the simplex, so the tolerance is on the parameter scale. A gradient
  rule is refused when the run starts.

- step:

  Size of the initial simplex, relative to each coordinate of the
  starting value. A single positive number, default `0.1`.

- adaptive:

  Use Gao and Han's dimension-dependent coefficients? `TRUE` or `FALSE`,
  default `TRUE`; see below. At \\p = 2\\ the two settings are
  identical.

- max_restarts:

  How many times a collapsed simplex may be rebuilt. A single
  non-negative number, default `3`; `0` turns the safeguard off.

- degenerate_tol:

  Rebuild the simplex when its conditioning falls below this, a single
  positive number, default `1e-6`. A non-positive value is refused, with
  a message naming `tol` rather than the argument.

- simplex:

  An optional starting simplex: a matrix with \\p+1\\ rows, one vertex
  per row. Defaults to `NULL`, meaning build one from `par` and `step`.
  Anything that is not a matrix or `NULL` is refused.

- maxit:

  Maximum iterations. Defaults to 2000.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`: no evaluation
  budget, so the run stops on the criterion or on `maxit`. Set a finite
  value to cap the cost of a run.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 50.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

## Value

An S7 object of class
[NelderMead](https://statmodels7.github.io/optimizers7/reference/NelderMead-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Details

Reflect the worst vertex through the centroid of the rest; if that is
the best point yet, try going further; if it is no better than the
second worst, pull back; and if even that fails, shrink everything
towards the best vertex. No model of the function is built and no
derivative is assumed to exist, which is why the method survives a kink.
It is also why it is slow: an ordering of \\p+1\\ values is very little
information about a surface, and on the quadratic below it costs 133
evaluations where
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
costs 3.

## Degenerate simplices

Nelder-Mead can converge to a point that is not a minimizer. McKinnon
(1998) exhibited a strictly convex function with continuous derivatives
on which it performs inside contractions for ever: the simplex flattens
onto a line through a point where the gradient is not zero, every vertex
agrees, and every ordinary stopping rule reports success. The *values*
offer no defense, behaving exactly as convergence would; what has gone
wrong is the *shape*.

So the shape is watched. The conditioning measured is \\\lvert \det E
\rvert\\ divided by the product of the edge lengths, where \\E\\ holds
the edges from the best vertex. It is 1 for a right-angled simplex, 0
for one that has collapsed into a lower dimension, and unchanged by
rescaling, so one threshold serves at every size. When it falls below
`degenerate_tol` the simplex is rebuilt right-angled at the current best
vertex, at the diameter it had reached, so the scale the run has earned
is kept. Each rebuild is counted, appears in the trace as `restart` and
is reported in the result's `message`.

The safeguard is not free: a restart costs \\p+1\\ evaluations and can
delay a genuine convergence. `max_restarts = 0` turns it off.

## Adaptive coefficients

The classical reflection, expansion, contraction and shrink factors are
\\1, 2, 1/2, 1/2\\, chosen when the method was proposed for two or three
parameters. In higher dimension a fixed expansion of 2 makes the simplex
overshoot along whichever direction it happened to try. Gao and Han
(2012) replace them by \\1,\\ 1 + 2/p,\\ 3/4 - 1/(2p),\\ 1 - 1/p\\,
which at \\p = 2\\ are *exactly* the classical values, so `TRUE` is the
default at no cost to the small problems anyone would recognize. At \\p
= 10\\ they are \\1, 1.2, 0.7, 0.9\\, and on a quadratic in ten unknowns
the adaptive run takes 959 iterations against 1230.

## When to reach for it

Rarely, and knowingly. On a smooth objective every gradient-based method
here beats it by orders of magnitude. Its place is an objective that is
genuinely non-smooth or noisy and whose subgradients are unavailable.
Where subgradients *are* available,
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
uses them and this method uses nothing.

## References

Nelder, J. A. and Mead, R. (1965). A simplex method for function
minimization. *The Computer Journal* **7**, 308–313.

McKinnon, K. I. M. (1998). Convergence of the Nelder-Mead simplex method
to a nonstationary point. *SIAM Journal on Optimization* **9**, 148–158.

Gao, F. and Han, L. (2012). Implementing the Nelder-Mead simplex
algorithm with adaptive parameters. *Computational Optimization and
Applications* **51**, 259–277.

## See also

[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
for the derivative-free method with a convergence theorem,
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
when subgradients are available,
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for the rule this method reads.

## Examples

``` r
nelder_mead()
#> <optimizer> nelder-mead
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 2000, evaluations Inf
#>   settings  : step = 0.1, adaptive = TRUE, max_restarts = 3, degenerate_tol = 1e-06, simplex = <NULL>

# A smooth problem, to show that it works and what it costs. BFGS reaches
# the same point in 3 evaluations.
q <- function(p) sum((p - c(1, 2))^2)
minimize(nelder_mead(), q, c(0, 0))@counts[["f"]]
#> [1] 133
minimize(bfgs(), q, c(0, 0), gr = function(p) 2 * (p - c(1, 2)))@counts[["f"]]
#> [1] 3

# What it is for: a sum of absolute deviations, whose minimizer is the
# median and whose derivative does not exist there.
set.seed(1)
y <- rcauchy(101)
r <- minimize(nelder_mead(), function(p) sum(abs(y - p)), par = 0)
c(nelder_mead = r@par, median = median(y))
#> nelder_mead      median 
#>  0.07342868  0.07342868 
abs(r@par - median(y))
#> [1] 2.477776e-10

# The adaptive coefficients are the classical ones at p = 2 and differ
# above it, where they save iterations.
f10 <- function(x) sum((x - 1:10)^2)
c(adaptive = minimize(nelder_mead(maxit = 20000), f10, numeric(10))@iterations,
  classical = minimize(nelder_mead(adaptive = FALSE, maxit = 20000),
                       f10, numeric(10))@iterations)
#>  adaptive classical 
#>       959      1230 

# The trace names the simplex operation taken at each iteration.
unique(minimize(nelder_mead(keep_trace = TRUE), q, c(0, 0))@trace$safeguard)
#> [1] "expand"       "reflect"      "contract out" "contract in" 
```
