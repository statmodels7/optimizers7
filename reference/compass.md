# Pattern Search, With Coordinate or Random Polling

Evaluates the objective around the current point along a set of
directions, moves to the first or to the best improvement, and shrinks
the radius when it finds none. It uses no derivative, and, unlike
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
it has a convergence theorem.

## Usage

``` r
compass(
  criterion = crit_stationary(),
  step = 0.1,
  directions = c("mads", "coordinate"),
  opportunistic = TRUE,
  expand = 2,
  shrink = 0.5,
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
  object. Defaults to `crit_stationary(1e-8)`, read on the poll size
  \\\Delta\\, the quantity in which the convergence theory is stated:
  for a continuously differentiable objective, the limit points of a
  coordinate pattern search with \\\Delta \to 0\\ are stationary
  (Torczon 1997), and for a locally Lipschitz objective the limit points
  of MADS are Clarke stationary (Audet and Dennis 2006).

- step:

  Initial poll size: `step` times the larger of 1 and the largest
  absolute coordinate of the starting value. A single positive number,
  default `0.1`.

- directions:

  `"mads"` (default) or `"coordinate"`. Partial matching applies, and
  any other string is rejected by
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

- opportunistic:

  Move to the first improvement found instead of polling every
  direction? `TRUE` or `FALSE`, default `TRUE`.

- expand:

  Factor applied to the poll size after a success. A single number of at
  least 1, default `2`.

- shrink:

  Factor applied after a failure. A single number strictly inside \\(0,
  1)\\, default `0.5`.

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
[Compass](https://statmodels7.github.io/optimizers7/reference/Compass-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Details

The directions form a *positive spanning set*: every vector in the space
is a non-negative combination of them. The theorem rests on this. If the
current point is not stationary, some direction in the set goes
downhill, so a poll that fails in every direction at a small radius
indicates that the point is close to stationary. A failed poll therefore
justifies shrinking the radius, and the limit points of a run whose
radius goes to zero are stationary.

## Poll directions

- `"coordinate"`:

  the \\2p\\ signed axes, which is compass search. This poll is cheap,
  deterministic and reproducible without a seed. Its convergence theorem
  assumes that \\f\\ is continuously differentiable.

- `"mads"`:

  a fresh random orthonormal basis at every poll, taken plus and minus.

The choice decides the outcome on the problems for which the method is
intended. Where \\f\\ is only Lipschitz, a *fixed* set of directions can
fail: if a kink has a diagonal ridge, no coordinate direction descends
it, so the poll fails at a point that is not stationary and the run
stops there. On \\\lvert x_1 + x_2 \rvert + 0.1\lVert x \rVert^2\\,
whose minimum is 0, coordinate polling from \\(1, 0.5)\\ stops at
`0.05`, whereas the random poll reaches values several orders of
magnitude lower, which depend on the seed (see the examples).

The remedy, which is the idea behind MADS, is to let the directions used
over the whole run become dense in the sphere, so that every direction
of descent is eventually polled; drawing a new orthonormal basis at each
poll achieves this with probability one. The implementation follows this
idea and is not LTMADS as published. The cost is reproducibility: a
random poll draws from the generator of R, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs the run and
the state is recorded in the `seed` of the result.

The poll size is capped at \\10^4\\ times its initial value, and the
trace records `poll capped` when the cap applies. A run whose poll size
falls below \\10^{-300}\\ stops with `converged = FALSE` and the message
`the poll size underflowed`.

## Opportunistic polling

Accepting the first improvement instead of the best gives a possibly
worse direction and saves up to \\2p - 1\\ evaluations per poll. The
default is `TRUE`; the last example compares the evaluation counts of
the two settings.

## References

Torczon, V. (1997). On the convergence of pattern search algorithms.
*SIAM Journal on Optimization* **7**, 1–25.

Audet, C. and Dennis, J. E. (2006). Mesh adaptive direct search
algorithms for constrained optimization. *SIAM Journal on Optimization*
**17**, 188–217.

## See also

[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
for the simplex,
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
when subgradients are available,
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for the default rule of this method.

## Examples

``` r
compass()
#> <optimizer> pattern search (mads)
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 2000, evaluations Inf
#>   settings  : step = 0.1, directions = mads, opportunistic = TRUE, expand = 2, shrink = 0.5

# A kink along a diagonal, which no coordinate direction descends. The
# true minimum is 0. Coordinate polling stops at the same wrong point every
# time; the random poll passes it.
f <- function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2)
minimize(compass(directions = "coordinate"), f, c(1, 0.5))@value
#> [1] 0.05
minimize(compass(directions = "coordinate"), f, c(1, 0.5))@value
#> [1] 0.05
set.seed(1)
minimize(compass(directions = "mads"), f, c(1, 0.5))@value
#> [1] 0.02685585

# A coordinate poll needs no seed; a random poll needs one, and records
# the state that it used.
a <- minimize(compass(directions = "coordinate"), f, c(1, 0.5))
b <- minimize(compass(directions = "coordinate"), f, c(1, 0.5))
c(identical(a@par, b@par), is.null(a@seed))
#> [1] TRUE TRUE

set.seed(7); m1 <- minimize(compass(), f, c(1, 0.5))
set.seed(7); m2 <- minimize(compass(), f, c(1, 0.5))
c(identical(m1@par, m2@par), length(m1@seed) > 0)
#> [1] TRUE TRUE

# Polling every direction costs more evaluations for the same answer.
q <- function(p) sum((p - c(1, 2))^2)
set.seed(5); first <- minimize(compass(), q, c(0, 0))
set.seed(5); every <- minimize(compass(opportunistic = FALSE), q, c(0, 0))
c(first = first@counts[["f"]], every = every@counts[["f"]])
#> first every 
#>   332   393 
```
