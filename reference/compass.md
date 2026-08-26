# Pattern Search, With Coordinate or Random Polling

Looks around the current point along a set of directions; moves to the
first or best improvement it finds, and shrinks the radius when it finds
none. Uses no derivative, and unlike
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
it comes with a convergence theorem.

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
  \\\Delta\\. This is the rule with a theorem behind it: the limit
  points of a pattern search with \\\Delta \to 0\\ are Clarke
  stationary.

- step:

  Initial poll size, scaled by the largest coordinate of the starting
  value. A single positive number, default `0.1`.

- directions:

  `"mads"` (default) or `"coordinate"`. Partial matching applies and any
  other string is refused, naming both.

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
is a non-negative combination of them. That is what buys the theorem. If
the current point is not stationary then some direction in the set goes
downhill, so a poll that fails everywhere is evidence about the point
and not about the directions. A failed poll therefore licenses shrinking
the radius, and the limit points of a run whose radius goes to zero are
stationary.

## Poll directions

- `"coordinate"`:

  the \\2p\\ signed axes, which is compass search. Cheap, deterministic
  and reproducible without a seed. The theorem it enjoys assumes \\f\\
  is continuously differentiable.

- `"mads"`:

  a fresh random orthonormal basis at every poll, taken plus and minus.

The difference decides the outcome on the problems this method exists
for. Where \\f\\ is merely Lipschitz, a *fixed* set of directions can
fail: a kink whose ridge runs diagonally is descended by no coordinate
direction, so the poll fails at a point that is not stationary and the
run stops there. Measured on \\\lvert x_1 + x_2 \rvert + 0.1\lVert x
\rVert^2\\, whose minimum is 0, coordinate polling from \\(1, 0.5)\\
stops at `0.05` on every one of ten runs, the same value each time; the
random poll has a median of `0.014` and reaches `9.8e-05` at best.

The repair, which is the idea behind MADS, is to let the directions used
over the whole run become dense in the sphere, so no direction of
descent is missed for ever, and drawing a new orthonormal basis at each
poll achieves that with probability one. What is implemented is that
idea rather than LTMADS as published, and it is the property the
convergence proof rests on. The cost is reproducibility: a random poll
draws from R's generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs the run and
the state is recorded in the result's `seed`.

## Opportunistic polling

Accepting the first improvement instead of the best costs a worse
direction and saves up to \\2p - 1\\ evaluations per poll. Measured on
the quadratic below: 332 evaluations against 393, in 96 iterations
against 98. The saving in evaluations is real and the cost in iterations
was not visible here, so `TRUE` is the default.

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
for the rule this method reads.

## Examples

``` r
compass()
#> <optimizer> pattern search (mads)
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 2000, evaluations Inf
#>   settings  : step = 0.1, directions = mads, opportunistic = TRUE, expand = 2, shrink = 0.5

# A kink running diagonally, which no coordinate direction descends. The
# true minimum is 0. Coordinate polling stops at the same wrong point every
# time; the random poll gets past it.
f <- function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2)
minimize(compass(directions = "coordinate"), f, c(1, 0.5))@value
#> [1] 0.05
minimize(compass(directions = "coordinate"), f, c(1, 0.5))@value
#> [1] 0.05
set.seed(1)
minimize(compass(directions = "mads"), f, c(1, 0.5))@value
#> [1] 0.02685585

# A coordinate poll needs no seed, a random one does, and it records the
# state it used.
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
