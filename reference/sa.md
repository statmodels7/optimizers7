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

  The stopping rule. Defaults to `crit_stationary(1e-8)`, read on the
  measure described above.

- visiting:

  `"uniform"` (default) or `"cauchy"`.

- t0:

  The initial temperature. `NULL`, the default, calibrates it from the
  objective.

- cooling:

  The factor the temperature is multiplied by at each level. Defaults to
  `0.85`.

- cycles:

  How many step-adjustment cycles per temperature level. Defaults to
  `3`.

- steps:

  How many sweeps of every coordinate per cycle. Defaults to `10`. One
  temperature level therefore costs `cycles * steps * length(par)`
  evaluations.

- step:

  The initial step, relative to the starting value. Defaults to `1`.

- target_accept:

  The acceptance rate the adaptation aims at, held inside a band of 0.1
  either side. Defaults to `0.5`.

- adjust:

  How hard the step is moved towards that rate. Defaults to `2`.

- n_eps:

  How many temperature levels the stationarity measure looks back over.
  Defaults to `4`.

- maxit:

  Maximum temperature levels. Defaults to 100. It is the budget that
  decides how tightly the run finishes: the step adaptation shrinks the
  proposal as the acceptance rate falls with the temperature, so on a
  quadratic the answer improves from \\2\times10^{-1}\\ at 15 levels to
  \\5\times10^{-5}\\ at 100.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many levels. Defaults to 10.

- keep_trace:

  Store the path? Defaults to `FALSE`.

## Value

An S7 object of class `Sa`, inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Details

The method exists here for the problems the local ones cannot start on:
a multimodal objective where the answer depends on which basin the run
began in. It is not a competitor to
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md) or
[`newton`](https://statmodels7.github.io/optimizers7/reference/newton.md)
on a smooth problem, where it will be beaten by orders of magnitude, and
the intended use is to hand its result to one of them — see
[`chain`](https://statmodels7.github.io/optimizers7/reference/chain.md).

**The adaptive step.** The parameters are moved one coordinate at a
time, and after every `steps` sweeps each coordinate's step is
multiplied or divided according to how often its moves were accepted, so
that the rate is held inside a band around `target_accept` (Corana et
al. 1987). This is what makes the method usable on a statistical
objective: the coordinates of an unconstrained parameter vector sit on
scales orders of magnitude apart, and one step length is wrong for all
of them. A coordinate accepting almost everything is being proposed too
timidly to explore; one accepting almost nothing is being thrown too far
to land.

**The proposal.** `"uniform"` draws the move uniformly on the
coordinate's current step, which with the adaptation above is Corana's
algorithm. `"cauchy"` draws it from a Cauchy, whose heavy tail lets a
run leave a basin in one move rather than walking out of it; that is
fast simulated annealing (Szu and Hartley 1987), and it is the \\q = 2\\
member of the Tsallis family. The general Tsallis visiting distribution
at arbitrary \\q\\ is deliberately not offered: its generator would have
to be transcribed and could not be checked against anything already
here, and an unverified generator is worse than a case that can be
verified.

**The temperature.** With `t0 = NULL` the initial temperature is
calibrated from the objective's own variation, by sampling proposals
around the starting value and setting \\T_0\\ so that an average uphill
move is accepted about four times in five. A fixed number cannot serve:
on an objective of size \\10^6\\ every move is accepted and the run is a
random walk, on one of size \\10^{-6}\\ none is and it is a poor local
search.

**What the run returns is the best point SEEN**, not the last one. An
annealing run wanders by construction, so its final iterate is a draw
and not an answer.

**Whether it converged is a separate question** and is never answered by
the schedule having finished. The stationarity measure reported is
Corana's own termination rule — by how much the best value has moved
over the last `n_eps` temperature levels — so
[`crit_stationary`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
IS that rule rather than a second convention invented beside it. A run
that merely exhausts `maxit` reports `converged = FALSE`, which for a
global search is the ordinary outcome and not a failure.

**The run is stochastic**, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) governs it and the
state it started from is recorded in the result.

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

[`chain`](https://statmodels7.github.io/optimizers7/reference/chain.md),
[`multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.md),
[`crit_stationary`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)

## Examples

``` r
sa()
#> <optimizer> simulated annealing (uniform)
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 100, evaluations Inf
#>   settings  : visiting = uniform, t0 = <NULL>, cooling = 0.85, cycles = 3, steps = 10, step = 1, target_accept = 0.5, adjust = 2, n_eps = 4

# an objective with many local minima, where a local method stops in the
# basin it started in and this one does not
rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
set.seed(1)
minimize(sa(), rastrigin, c(4.4, -3.6))@value
#> [1] 0.9141033
minimize(bfgs(), rastrigin, c(4.4, -3.6))@value
#> [1] 33.82832
```
