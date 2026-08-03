# Pattern Search, With Coordinate or Random Polling

Looks around the current point along a set of directions; moves to the
first or best improvement it finds, and shrinks the radius when it finds
none. Uses no derivative, and unlike
[`nelder_mead`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
it comes with a convergence theorem.

## Usage

``` r
compass(
  criterion = crit_stationary(1e-08),
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

  The stopping rule. Defaults to `crit_stationary(1e-8)`, on the poll
  size.

- step:

  Initial poll size, scaled by the largest coordinate of the starting
  value. Defaults to `0.1`.

- directions:

  `"mads"` (default) or `"coordinate"`; see Details.

- opportunistic:

  Move to the first improvement found rather than polling every
  direction? Defaults to `TRUE`.

- expand:

  Factor applied to the poll size after a success. Defaults to `2`.

- shrink:

  Factor applied after a failure. Defaults to `0.5`.

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

An S7 object of class `Compass`, inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Details

The directions form a *positive spanning set*: any vector in the space
is a non-negative combination of them. That is the whole idea, and it is
what buys the theorem — if the current point is not stationary then some
direction in the set goes downhill, so a poll that fails everywhere is
evidence about the point rather than about the directions. A failed poll
therefore licenses shrinking the radius, and the limit points of a run
in which the radius goes to zero are stationary.

### Poll directions

- `"coordinate"`:

  the \\2p\\ signed axes: this is compass search. Cheap, deterministic,
  reproducible. The theorem it enjoys assumes \\f\\ is continuously
  differentiable.

- `"mads"`:

  a fresh random orthonormal basis at every poll, taken plus and minus.

The difference is not cosmetic on the problems this method exists for.
When \\f\\ is merely Lipschitz — has kinks — a *fixed* set of directions
can fail: a kink whose ridge runs diagonally is descended by no
coordinate direction, the poll fails at a point that is not stationary,
and the run stops there. The repair, which is the idea behind MADS, is
to let the set of directions used over the whole run become dense in the
sphere, so no direction of descent can be missed for ever; drawing a new
orthonormal basis at each poll achieves that with probability one. What
is implemented here is that idea rather than LTMADS as published, but it
is the property the convergence proof rests on.

The cost is reproducibility: a random poll draws from R's generator, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs the run.

### Opportunistic polling

Accepting the first improvement rather than the best costs a worse
direction and saves up to \\2p - 1\\ evaluations. On an expensive
objective that trade is usually favourable, so it is the default; on a
cheap one, `opportunistic = FALSE` tends to need fewer iterations.

## References

Torczon, V. (1997). On the convergence of pattern search algorithms.
*SIAM Journal on Optimization* **7**, 1–25.

Audet, C. and Dennis, J. E. (2006). Mesh adaptive direct search
algorithms for constrained optimization. *SIAM Journal on Optimization*
**17**, 188–217.

## See also

[`nelder_mead`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
[`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md),
[`crit_stationary`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)

## Examples

``` r
compass()
#> <optimizer> pattern search (mads)
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 2000, evaluations Inf
#>   settings  : step = 0.1, directions = mads, opportunistic = TRUE, expand = 2, shrink = 0.5

# a kink running diagonally, which no coordinate direction descends
f <- function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2)
minimize(compass(), f, c(1, 0.5))@value
#> [1] 0.002114785
set.seed(1)
minimize(compass(directions = "mads"), f, c(1, 0.5))@value
#> [1] 0.02685585
```
