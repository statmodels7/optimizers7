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
  object. Defaults to `crit_stationary(1e-8)`, read on the size of the
  simplex (the largest distance from the best vertex to another vertex),
  so the tolerance is on the parameter scale. A gradient rule is
  rejected when the run starts.

- step:

  Size of the initial simplex: the vertex that moves coordinate \\j\\ is
  displaced by `step` times the larger of 1 and the absolute value of
  that coordinate of the starting value. A single positive number,
  default `0.1`.

- adaptive:

  Use Gao and Han's dimension-dependent coefficients? `TRUE` or `FALSE`,
  default `TRUE`; see below. At \\p = 2\\ the two settings are
  identical.

- max_restarts:

  How many times a collapsed simplex may be rebuilt. A single
  non-negative number, default `3`; `0` turns the safeguard off.

- degenerate_tol:

  Rebuild the simplex when its conditioning falls below this, a single
  positive number, default `1e-6`. A value that is not positive is
  rejected.

- simplex:

  An optional starting simplex: a matrix with \\p+1\\ rows, one vertex
  per row. Defaults to `NULL`, meaning build one from `par` and `step`.
  A value that is not a matrix or `NULL` is rejected.

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
derivative is assumed to exist, so the method works at a kink. It is
slow because the ordering of \\p+1\\ values carries little information
about the surface; the first example compares its evaluation count with
that of
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md).

## Degenerate simplices

Nelder-Mead can converge to a point that is not a minimizer. McKinnon
(1998) exhibited a strictly convex function with continuous derivatives
on which it performs inside contractions indefinitely: the simplex
flattens onto a line through a point where the gradient is not zero, the
vertices approach a common value, and every ordinary stopping rule
reports success. The function values behave as they would at
convergence, and only the shape of the simplex shows the failure.

The method therefore monitors the shape. The conditioning measured is
\\\lvert \det E \rvert\\ divided by the product of the edge lengths,
where \\E\\ holds the edges from the best vertex. It is 1 for a
right-angled simplex, 0 for one that has collapsed into a lower
dimension, and unchanged by rescaling, so one threshold serves at every
size. When it falls below `degenerate_tol` the simplex is rebuilt with
right angles at the current best vertex, with edges equal to the size
that the simplex had reached. Each rebuild is counted, appears in the
trace as `restart` and is reported in the `message` of the result.

The safeguard is not free: a restart costs \\p+1\\ evaluations and can
delay a genuine convergence. `max_restarts = 0` turns it off.

## Adaptive coefficients

The classical reflection, expansion, contraction and shrink factors are
\\1, 2, 1/2, 1/2\\, chosen when the method was proposed for two or three
parameters. In higher dimension a fixed expansion of 2 makes the simplex
overshoot along whichever direction it happened to try. Gao and Han
(2012) replace them by \\1,\\ 1 + 2/p,\\ 3/4 - 1/(2p),\\ 1 - 1/p\\,
which at \\p = 2\\ are *exactly* the classical values, so `TRUE` is the
default without changing the method on problems with two parameters. At
\\p = 10\\ they are \\1, 1.2, 0.7, 0.9\\, and on a quadratic in ten
unknowns the adaptive run needs fewer iterations than the classical one
(see the examples).

## Where the method applies

On a smooth objective the gradient-based methods that use curvature,
such as
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md),
need far fewer evaluations. The method is intended for an objective that
is non-smooth or noisy and whose subgradients are unavailable. Where
subgradients *are* available,
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
uses them.

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
for the default rule of this method.

## Examples

``` r
nelder_mead()
#> <optimizer> nelder-mead
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 2000, evaluations Inf
#>   settings  : step = 0.1, adaptive = TRUE, max_restarts = 3, degenerate_tol = 1e-06, simplex = <NULL>

# A smooth problem, to show the cost: BFGS reaches the same point in far
# fewer evaluations.
q <- function(p) sum((p - c(1, 2))^2)
minimize(nelder_mead(), q, c(0, 0))@counts[["f"]]
#> [1] 133
minimize(bfgs(), q, c(0, 0), gr = function(p) 2 * (p - c(1, 2)))@counts[["f"]]
#> [1] 3

# A sum of absolute deviations, whose minimizer is the median and whose
# derivative does not exist there.
set.seed(1)
y <- rcauchy(101)
r <- minimize(nelder_mead(), function(p) sum(abs(y - p)), par = 0)
c(nelder_mead = r@par, median = median(y))
#> nelder_mead      median 
#>  0.07342868  0.07342868 
abs(r@par - median(y))
#> [1] 2.477776e-10

# The adaptive coefficients equal the classical ones at p = 2 and differ
# above it; here they save iterations.
f10 <- function(x) sum((x - 1:10)^2)
c(adaptive = minimize(nelder_mead(maxit = 20000), f10, numeric(10))@iterations,
  classical = minimize(nelder_mead(adaptive = FALSE, maxit = 20000),
                       f10, numeric(10))@iterations)
#>  adaptive classical 
#>       969      1230 

# The trace names the simplex operation taken at each iteration.
unique(minimize(nelder_mead(keep_trace = TRUE), q, c(0, 0))@trace$safeguard)
#> [1] "expand"       "reflect"      "contract out" "contract in" 
```
