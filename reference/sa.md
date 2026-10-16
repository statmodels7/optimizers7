# Simulated Annealing With an Adaptive Step

Creates an optimizer for a global search that accepts uphill moves with
a probability that falls as the run cools, with one step length per
coordinate that is adjusted during the run to hold the acceptance rate
of that coordinate near a target.

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
  above. A rule that reads a gradient is rejected when the run starts,
  because this method computes none.

- visiting:

  `"uniform"` (default) or `"cauchy"`. Partial matching applies, and any
  other string is rejected by
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

- t0:

  The initial temperature, a single positive number, or `NULL` (the
  default) to calibrate it from the objective at a cost of
  `max(10, 10 * length(par))` evaluations.

- cooling:

  The factor by which the temperature is multiplied at each level, a
  single number strictly inside \\(0, 1)\\. Defaults to `0.85`. The
  value `1` would never cool and `0` would freeze at once, so both
  endpoints are rejected.

- cycles:

  How many step-adjustment cycles per temperature level, a positive
  whole number. Defaults to `3`.

- steps:

  How many sweeps of every coordinate per cycle, a positive whole
  number. Defaults to `10`. One temperature level therefore costs
  `cycles * steps * length(par)` evaluations, and a whole run costs that
  times `maxit`, plus one for the starting value and
  `max(10, 10 * length(par))` more if the temperature is calibrated.

- step:

  The initial step: the step of coordinate \\j\\ is `step` times the
  larger of 1 and the absolute value of that coordinate of the starting
  value. Defaults to `1`. The adaptation changes the step of each
  coordinate within the first few cycles.

- target_accept:

  The target acceptance rate of the adaptation, held inside a band of
  0.1 on either side. A single number strictly inside \\(0.1, 0.9)\\; a
  value outside is rejected, because the band around it would reach 0
  or 1. Defaults to `0.5`.

- adjust:

  How strongly the step is moved towards that rate, a single positive
  number. Defaults to `2`, which is Corana's constant.

- n_eps:

  How many temperature levels the stationarity measure covers, a
  positive whole number. Defaults to `4`.

- maxit:

  Maximum temperature levels. Defaults to 100. The number of levels sets
  how close to the minimizer the run finishes, because the step
  adaptation shrinks the proposal as the acceptance rate falls with the
  temperature.

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

The method is intended for problems on which a local method depends on
its start: a multimodal objective, on which the answer depends on the
basin in which the run begins. On a smooth problem
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
need far fewer evaluations, and the intended use is to pass the result
of this method to one of them; see
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md).

## The adaptive step

The parameters are moved one coordinate at a time, and after every
`steps` sweeps the step of each coordinate is multiplied or divided
according to how often its moves were accepted, so that the rate is held
inside a band around `target_accept` (Corana et al. 1987). The
adaptation makes the method usable on a statistical objective, whose
unconstrained coordinates lie on scales orders of magnitude apart, so
that no single step length suits all of them. A coordinate whose moves
are almost all accepted has a step that is too short to explore, and a
coordinate whose moves are almost all rejected has a step that is too
long.

## The proposal

`"uniform"` draws the move uniformly on the current step of the
coordinate, which together with the adaptation above is Corana's
algorithm. `"cauchy"` draws it from a Cauchy distribution, whose heavy
tail lets a run leave a basin in one move. This is the visiting
distribution of fast simulated annealing (Szu and Hartley 1987), the \\q
= 2\\ member of the Tsallis family; the temperature still cools
geometrically, as for the uniform proposal.

On a surface with many basins the Cauchy proposal reaches low values
more often than the uniform one; the examples compare the two on the
Rastrigin function. The default is nevertheless `"uniform"`, which is
Corana's algorithm as published and the better-behaved proposal where
the objective is defined on a bounded region.

The general Tsallis visiting distribution at arbitrary \\q\\ is not
implemented; only the uniform and the Cauchy proposals are.

## The temperature

With `t0 = NULL` the initial temperature is calibrated from the
variation of the objective: the method samples `max(10, 10 * p)`
proposals around the starting value, where \\p\\ is the number of
parameters, averages the absolute change of the objective over them, and
sets \\T_0\\ so that a change of that average size is accepted with
probability 0.8. The calibration costs those evaluations once, before
the run begins.

A fixed temperature does not suit every objective, because the
Metropolis probability \\\exp(-\Delta f / T)\\ compares the temperature
with the scale of the objective. At \\T = 1\\ a step that raises the
objective by \\10^{-6}\\ is accepted with probability close to 1 and the
run is a random walk; one that raises it by \\10^{6}\\ is accepted with
probability close to 0, and the method is a local search. For an
objective of small scale, such as \\10^{-6}\sum(x - (1,2))^2\\, a
calibrated temperature therefore gets closer to the minimizer than
`t0 = 1`.

On a convex objective a low temperature is not a handicap, since there
is a single basin, and the calibration can leave the run further from
the minimizer than a fixed `t0 = 1` within the same budget.

## The value returned and the meaning of convergence

The result is the best point seen, not the last one, since the last
point of an annealing run is a random draw from the walk. The `value`
column of the trace is correspondingly the best value so far, and is
monotone non-increasing.

The end of the schedule does not mean convergence. The stationarity
measure is Corana's termination rule. At the end of each temperature
level it is the largest absolute difference between the value at which
the walk stopped and each of two references: the best value found so
far, and the values at which the walk stopped at the last `n_eps`
levels. It is infinite until `n_eps + 1` levels have run, and
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
reads it. A run that exhausts `maxit` reports `converged = FALSE`, which
for a global search is the ordinary outcome.

The run is stochastic, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs it, and the
state from which it began is recorded in the `seed` of the result, so
that it can be repeated.

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
for the rule that reads the measure of this method.

## Examples

``` r
sa()
#> <optimizer> simulated annealing (uniform)
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 100, evaluations Inf
#>   settings  : visiting = uniform, t0 = <NULL>, cooling = 0.85, cycles = 3, steps = 10, step = 1, target_accept = 0.5, adjust = 2, n_eps = 4

# Rastrigin has many local minima (121 in the box [-5.12, 5.12]^2). A
# local method stops in the basin in which it started; this one leaves it.
rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
set.seed(1)
minimize(sa(), rastrigin, c(4.4, -3.6))@value    # 0.914
#> [1] 0.9141033
minimize(bfgs(), rastrigin, c(4.4, -3.6))@value  # 7.96
#> [1] 7.959662

# The heavy-tailed proposal leaves a basin in one move, and on this surface
# that gives a value about two orders of magnitude lower.
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

# A gradient rule is rejected when the run starts, since this method
# computes no gradient.
try(minimize(sa(criterion = crit_grad()), rastrigin, c(4.4, -3.6)))
#> Error : The stopping rule needs gradient, which simulated annealing (uniform) does not provide.
#>   Choose a criterion that this optimizer can evaluate, or a method that provides it.
```
