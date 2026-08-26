# Simulated Annealing With an Adaptive Step

A global search that accepts uphill moves with a probability falling as
the run cools, with one step length per coordinate adjusted in flight to
hold that coordinate's acceptance rate near a target.

## Usage

``` r
sa(
  criterion = crit_stationary(),
  visiting = c("uniform", "cauchy"),
  t0 = NULL,
  cooling = 0.85,
  cycles = 3,
  steps = 10,
  step = 1,
  target_accept = 0.5,
  adjust = 2,
  n_eps = 4,
  maxit = 100,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 10,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object. Defaults to `crit_stationary(1e-8)`, read on Corana's measure
  above. A rule reading a gradient is refused when the run starts, this
  method computing none.

- visiting:

  `"uniform"` (default) or `"cauchy"`. Partial matching applies, and any
  other string is refused by name.

- t0:

  The initial temperature, a single positive number, or `NULL` (the
  default) to calibrate it from the objective at a cost of 20
  evaluations.

- cooling:

  The factor the temperature is multiplied by at each level, a single
  number **strictly** inside \\(0, 1)\\. Defaults to `0.85`. `1` would
  never cool and `0` would freeze at once, so both endpoints are
  refused.

- cycles:

  How many step-adjustment cycles per temperature level, a positive
  whole number. Defaults to `3`.

- steps:

  How many sweeps of every coordinate per cycle, a positive whole
  number. Defaults to `10`. One temperature level therefore costs
  `cycles * steps * length(par)` evaluations, and a whole run costs that
  times `maxit`, plus one for the starting value and 20 more if the
  temperature is calibrated.

- step:

  The initial step, relative to the starting value. Defaults to `1`. It
  is only a starting point: the adaptation moves each coordinate's step
  from here within the first few cycles.

- target_accept:

  The acceptance rate the adaptation aims at, held inside a band of 0.1
  either side. A single number strictly inside \\(0.1, 0.9)\\; anything
  outside is refused, an extreme target being unreachable by any step.
  Defaults to `0.5`.

- adjust:

  How hard the step is moved towards that rate, a single positive
  number. Defaults to `2`, which is Corana's own constant.

- n_eps:

  How many temperature levels the stationarity measure looks back over,
  a positive whole number. Defaults to `4`.

- maxit:

  Maximum temperature levels. Defaults to 100. This is the budget that
  decides how tightly the run finishes, the step adaptation shrinking
  the proposal as the acceptance rate falls with the temperature.
  Measured on \\\sum(x - (1, 2))^2\\ from the origin at one seed, the
  distance to the answer goes `2.7e-02`, `1.8e-02`, `3.0e-03`, `2.3e-05`
  at 15, 30, 60 and 100 levels, for 921, 1821, 3621 and 6021
  evaluations.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many levels. Defaults to 10.

- keep_trace:

  Store the path? Defaults to `FALSE`. The `value` column is the best so
  far.

## Value

An S7 object of class
[Sa](https://statmodels7.github.io/optimizers7/reference/Sa-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Details

The method exists here for the problems the local ones cannot start on:
a multimodal objective where the answer depends on which basin the run
began in. It is not a competitor to
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
or
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
on a smooth problem, where it will be beaten by orders of magnitude, and
the intended use is to hand its result to one of them; see
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md).

## The adaptive step

The parameters are moved one coordinate at a time, and after every
`steps` sweeps each coordinate's step is multiplied or divided according
to how often its moves were accepted, so that the rate is held inside a
band around `target_accept` (Corana et al. 1987). The adaptation makes
the method usable on a statistical objective, whose unconstrained
coordinates sit on scales orders of magnitude apart: one step length is
wrong for all of them. A coordinate accepting almost everything is being
proposed too timidly to explore, and one accepting almost nothing is
being thrown too far to land.

## The proposal

`"uniform"` draws the move uniformly on the coordinate's current step,
which together with the adaptation above is Corana's algorithm.
`"cauchy"` draws it from a Cauchy, whose heavy tail lets a run leave a
basin in one move instead of walking out of it. That is fast simulated
annealing (Szu and Hartley 1987), the \\q = 2\\ member of the Tsallis
family.

The difference is large on a surface with many basins. Over thirty seeds
on Rastrigin in two dimensions from \\(4.4, -3.6)\\, with 40 temperature
levels: the Cauchy proposal has a median value of `0.0051` and gets
below `0.01` on 19 of 30 runs, the uniform one a median of `0.94` and 1
of 30. The default is nevertheless `"uniform"`, that being Corana's
algorithm as published and the better-behaved proposal where the
objective is defined on a bounded region.

The general Tsallis visiting distribution at arbitrary \\q\\ is not
offered. Its generator would have to be transcribed and could not be
checked against anything already here.

## The temperature

With `t0 = NULL` the initial temperature is calibrated from the
objective's own variation, by sampling proposals around the starting
value and setting \\T_0\\ so that an average uphill move is accepted
about four times in five. That costs 20 evaluations, once, before the
run begins.

A fixed number cannot serve, because the Metropolis probability
\\\exp(-\Delta f / T)\\ compares the temperature against the objective's
own scale. At \\T = 1\\ a step costing \\10^{-6}\\ is accepted with
probability 1 and the run is a random walk; one costing \\10^{6}\\ is
accepted with probability 0 and the method has become a local search.
Measured over eight seeds on \\10^{-6}\sum(x - (1,2))^2\\, the
calibrated run ends a median `0.0044` from the answer and `t0 = 1` a
median `0.375`.

The other side of that is worth knowing too: on a **convex** objective
too cold a temperature costs nothing, there being one basin to find. On
\\10^{6}\sum(x - (1,2))^2\\ the fixed `t0 = 1` reaches `4.6e-07` where
the calibration reaches `0.0013`. The calibration is insurance against a
multimodal surface, and on a problem that has no second basin it is a
small tax.

## What the run returns, and what convergence means here

The result is the best point **seen**, not the last one. An annealing
run wanders by construction, so its final iterate is a draw. The `value`
column of the trace is correspondingly the best so far, and is monotone
non-increasing.

Whether the run converged is a separate question, never answered by the
schedule having finished. The stationarity measure reported is Corana's
own termination rule, by how much the best value has moved over the last
`n_eps` temperature levels, so
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
is that rule and no second convention invented beside it. A run that
merely exhausts `maxit` reports `converged = FALSE`, which for a global
search is the ordinary outcome.

The run is stochastic, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs it, and the
state it began from is recorded in the result's `seed` for repeating it.

## References

Corana, A., Marchesi, M., Martini, C. and Ridella, S. (1987). Minimizing
multimodal functions of continuous variables with the simulated
annealing algorithm. *ACM Transactions on Mathematical Software* **13**,
262–280.

Kirkpatrick, S., Gelatt, C. D. and Vecchi, M. P. (1983). Optimization by
simulated annealing. *Science* **220**, 671–680.

Szu, H. and Hartley, R. (1987). Fast simulated annealing. *Physics
Letters A* **122**, 157–162.

## See also

[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
for handing the answer to a local method,
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
for the other way of covering several basins,
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for the rule this method's measure feeds.

## Examples

``` r
sa()
#> <optimizer> simulated annealing (uniform)
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 100, evaluations Inf
#>   settings  : visiting = uniform, t0 = <NULL>, cooling = 0.85, cycles = 3, steps = 10, step = 1, target_accept = 0.5, adjust = 2, n_eps = 4

# Rastrigin has 121 local minima. A local method stops in the basin it
# started in; this one leaves it.
rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
set.seed(1)
minimize(sa(), rastrigin, c(4.4, -3.6))@value    # 0.914
#> [1] 0.9141033
minimize(bfgs(), rastrigin, c(4.4, -3.6))@value  # 7.96
#> [1] 7.959662

# The heavy-tailed proposal leaves a basin in one move, and on this surface
# that is worth two orders of magnitude.
set.seed(2)
minimize(sa(visiting = "cauchy", maxit = 40), rastrigin, c(4.4, -3.6))@value
#> [1] 0.002005912
set.seed(2)
minimize(sa(visiting = "uniform", maxit = 40), rastrigin, c(4.4, -3.6))@value
#> [1] 0.5458509

# The reported value is the best seen, and the trace records the same:
# it never rises, though the walk itself does.
set.seed(5)
r <- minimize(sa(maxit = 30, keep_trace = TRUE), rastrigin, c(4.4, -3.6))
all(diff(r@trace$value) <= 0)
#> [1] TRUE
all.equal(r@value, min(r@trace$value))
#> [1] TRUE

# The run repeats exactly from the state it recorded.
set.seed(11)
a <- minimize(sa(maxit = 20), rastrigin, c(4.4, -3.6))
assign(".Random.seed", a@seed, globalenv())
b <- minimize(sa(maxit = 20), rastrigin, c(4.4, -3.6))
identical(a@par, b@par)
#> [1] TRUE

# A gradient rule is refused when the run starts, this method computing none.
try(minimize(sa(criterion = crit_grad()), rastrigin, c(4.4, -3.6)))
#> Error : The stopping rule needs gradient, which simulated annealing (uniform) does not provide.
#>   Choose a criterion this optimizer can evaluate, or a method that provides it.
```
