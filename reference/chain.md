# Run Optimizers One After Another

Builds one optimizer out of several, run in order, each starting from
the point the previous one reached. The composition a global search
needs: `chain(sa(), lbfgs())` explores and then descends from wherever
the exploration left off, and neither method has to know about the
other. The result is itself an
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
so it goes anywhere a single method goes,
`multistart(chain(sa(), lbfgs()))` included.

## Usage

``` r
chain(..., verbose = FALSE, keep_trace = FALSE)
```

## Arguments

- ...:

  Optimizers, in the order they should run. At least one is required;
  one is accepted. Anything that is not an
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  raises an error naming
  [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md)
  and
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  as examples.

- verbose:

  `TRUE` or `FALSE`, default `FALSE`: report which stage is starting,
  one line each. Independent of the stages' own `verbose`.

- keep_trace:

  `TRUE` or `FALSE`, default `FALSE`: assemble a trace across the
  stages. The stages must also have been built with `keep_trace = TRUE`
  to contribute anything.

## Value

An S7 object of class
[Chain](https://statmodels7.github.io/optimizers7/reference/Chain-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).
Its `name` is the stages' names joined by `then`, and its `criterion`,
`maxit` and `max_eval` are copied from the last stage.

## Each stage keeps its own settings

A stage carries its own stopping rule and its own budgets, which is the
reason a chain is built out of optimizers instead of taking arguments: a
coarse rule and a small budget for the exploration, a tight rule for the
descent, written as
`chain(sa(maxit = 200), newton(criterion = crit_grad(1e-12)))`.

## What the result reports

The point, the value, the gradient and `iterations` are the last
stage's, that being where the run ended. So is `converged`: a chain has
converged when the method that finished it says so. An earlier stage
exhausting its own budget is the ordinary way a global search ends and
leaves the flag alone.

The evaluation counts are **summed over the stages**, so a chain can be
compared against a single optimizer on total work. `message` is the
stages' messages joined, each prefixed with its stage number.

A trace is assembled when both the chain and the stage were built with
`keep_trace = TRUE`; the chain's flag decides whether one is kept at all
and each stage's decides whether that stage contributes rows. It carries
a `stage` column. Stages need not report the same columns, a
derivative-free method having a `stationarity` where a descent has a
`gnorm`; when they disagree only the last stage's trace is kept.

A stage that raises propagates: a method that cannot run on the
objective is a fact about the objective. A stage that runs without
converging simply passes its point on, and the first stage of a chain
usually does exactly that. Every stage's rule is checked by
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
before any of them runs, so a chain whose third stage cannot evaluate
its criterion fails before spending the first two.

`chain(x)` with a single stage is that stage's run reported through the
chain, identical to the stage's own down to the evaluation counts.

## See also

[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
for the other wrapper,
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) for
the search a chain usually opens with,
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
for the rule check.

## Examples

``` r
chain(sa(maxit = 10), bfgs())
#> <optimizer> simulated annealing (uniform) then BFGS
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 500, evaluations Inf
#>   settings  : stages = <list>

# Rastrigin has 121 local minima, so which one a descent method finds is
# decided by where it began. The search picks the basin; the descent
# finishes the job inside it.
rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
rg <- function(p) 2 * p + 20 * pi * sin(2 * pi * p)

set.seed(3)
minimize(bfgs(), rastrigin, c(3.5, -2.5), gr = rg)@value   # 12.93
#> [1] 12.93443
set.seed(3)
minimize(sa(maxit = 20), rastrigin, c(3.5, -2.5))@value    # 1.48
#> [1] 1.479952
set.seed(3)
minimize(chain(sa(maxit = 20), bfgs()), rastrigin, c(3.5, -2.5))@value
#> [1] 0.9949591

# The counts are the whole chain's work; the iteration count is the last
# stage's.
set.seed(3)
res <- minimize(chain(sa(maxit = 20), bfgs()), rastrigin, c(3.5, -2.5))
unlist(res@counts)
#>    f    g    h 
#> 1250    0    0 
res@iterations
#> [1] 4

# One stage is not a special case: it reports exactly what the stage does.
identical(minimize(chain(bfgs()), rastrigin, c(3.5, -2.5), gr = rg)@par,
          minimize(bfgs(), rastrigin, c(3.5, -2.5), gr = rg)@par)
#> [1] TRUE

# A stage that cannot evaluate its own rule is caught before anything runs.
try(minimize(chain(bfgs(), nelder_mead(criterion = crit_grad())),
             rastrigin, c(1, 1), gr = rg))
#> Error : The stopping rule needs gradient, which nelder-mead does not provide.
#>   Choose a criterion this optimizer can evaluate, or a method that provides it.
```
