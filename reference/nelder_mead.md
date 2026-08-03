# The Nelder-Mead Simplex Method

Keeps \\p+1\\ points, and at each step reflects the worst of them
through the face of the others. It uses no derivative and asks nothing
of the objective but that it return a number.

## Usage

``` r
nelder_mead(
  criterion = crit_stationary(1e-08),
  step = 0.1,
  adaptive = TRUE,
  max_restarts = 3,
  degenerate_tol = 1e-06,
  simplex = NULL,
  maxit = 2000,
  max_eval = 20000,
  verbose = FALSE,
  refresh = 50,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule. Defaults to `crit_stationary(1e-8)`, on the
  diameter of the simplex.

- step:

  Size of the initial simplex, relative to each coordinate of the
  starting value. Defaults to `0.1`.

- adaptive:

  Use Gao and Han's dimension-dependent coefficients? Defaults to
  `TRUE`; see Details.

- max_restarts:

  How many times a collapsed simplex may be rebuilt. Defaults to `3`.
  Zero disables the safeguard.

- degenerate_tol:

  Rebuild the simplex when its conditioning falls below this. Defaults
  to `1e-6`.

- simplex:

  An optional starting simplex: a matrix with \\p+1\\ rows, one vertex
  per row. Defaults to `NULL`, meaning build one from `par` and `step`.

- maxit:

  Maximum iterations. Defaults to 2000.

- max_eval:

  Maximum objective evaluations. Defaults to 20000.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 50.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

## Value

An S7 object of class `NelderMead`, inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Details

Reflect the worst vertex through the centroid of the rest; if that is
the best point yet, try going further; if it is no better than the
second worst, pull back; and if even that fails, shrink everything
towards the best vertex. No model of the function is built and no
derivative is assumed to exist, which is why it survives a kink — and
also why it is slow, since an ordering of \\p+1\\ values is very little
information about a surface.

### Degenerate simplices

Nelder-Mead can converge to a point that is not a minimiser. McKinnon
(1998) exhibited a strictly convex function with continuous derivatives
on which it performs inside contractions for ever: the simplex flattens
onto a line through a point where the gradient is not zero, every vertex
agrees, and every ordinary stopping rule reports success. There is no
defence in the *values* — they behave exactly as convergence would —
because what has gone wrong is the *shape* of the simplex.

So the shape is what is watched. The conditioning measured is \\\lvert
\det E \rvert\\ divided by the product of the edge lengths, where \\E\\
holds the edges from the best vertex: it is 1 for a right-angled
simplex, 0 for one that has collapsed into a lower dimension, and
unchanged by rescaling, so one threshold serves at every size. When it
falls below `degenerate_tol` the simplex is rebuilt right-angled at the
current best vertex, at the diameter it had reached — keeping the scale
the run has earned rather than starting over. Each rebuild is counted,
reported in the trace as `"restart"` and in the result's message.

The safeguard is not free: a restart costs \\p+1\\ evaluations and can
delay an honest convergence. `max_restarts = 0` turns it off.

### Adaptive coefficients

The classical reflection, expansion, contraction and shrink factors are
\\1, 2, 1/2, 1/2\\, chosen when the method was proposed for two or three
parameters. In higher dimension a fixed expansion of 2 makes the simplex
overshoot along whichever direction it happened to try. Gao and Han
(2012) replace them by \\1,\\ 1 + 2/p,\\ 3/4 - 1/(2p),\\ 1 - 1/p\\,
which at \\p = 2\\ reduce *exactly* to the classical values — so the
default is `TRUE` at no cost to the small problems anyone would
recognise.

### Scope

Rarely, and knowingly. If the objective is smooth, every gradient-based
method here will beat it by orders of magnitude. Its place is an
objective that is genuinely non-smooth or noisy and whose subgradients
are not available; if they *are* available,
[`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
is the better tool, since it uses them and this uses nothing.

## References

Nelder, J. A. and Mead, R. (1965). A simplex method for function
minimization. *The Computer Journal* **7**, 308–313.

McKinnon, K. I. M. (1998). Convergence of the Nelder-Mead simplex method
to a nonstationary point. *SIAM Journal on Optimization* **9**, 148–158.

Gao, F. and Han, L. (2012). Implementing the Nelder-Mead simplex
algorithm with adaptive parameters. *Computational Optimization and
Applications* **51**, 259–277.

## See also

[`compass`](https://statmodels7.github.io/optimizers7/reference/compass.md),
[`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md),
[`crit_stationary`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)

## Examples

``` r
nelder_mead()
#> <optimizer> nelder-mead
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 2000, evaluations 20000
#>   settings  : step = 0.1, adaptive = TRUE, max_restarts = 3, degenerate_tol = 1e-06, simplex = <NULL>

# a smooth problem, to show it works and to show what it costs
minimize(nelder_mead(), function(p) sum((p - c(1, 2))^2), c(0, 0))
#> <optimizer_result> nelder-mead
#>   value      : 1.66933e-18
#>   par        : 1 2
#>   iterations : 68   evaluations: f 133, g 0
#>   converged  : yes (stationarity < 1e-08)

# what it is actually for: a sum of absolute deviations, whose minimiser is
# the median and whose derivative does not exist there
set.seed(1)
y <- rcauchy(101)
minimize(nelder_mead(), function(p) sum(abs(y - p)), par = 0)@par
#> [1] 0.07342868
median(y)
#> [1] 0.07342868
```
