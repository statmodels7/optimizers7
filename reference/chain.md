# Run Optimizers One After Another

Builds one optimizer out of several optimizers that run in order, and
each of them starts from the point that the previous one reached. This
is the composition used with a global search: `chain(sa(), lbfgs())`
explores and then descends from the point where the exploration stopped,
and each method works without code written for the other. The result is
itself an
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
so it can be used wherever a single method can, including
`multistart(chain(sa(), lbfgs()))`.

## Usage

``` r
chain(..., verbose = FALSE, keep_trace = FALSE)
```

## Arguments

- ...:

  Optimizers, in the order in which they should run. At least one is
  required; one is accepted. Anything that is not an
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  raises an error naming
  [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md)
  and
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  as examples.

- verbose:

  `TRUE` or `FALSE`, default `FALSE`: report which stage is starting,
  one line each. It does not depend on the `verbose` of the stages.

- keep_trace:

  `TRUE` or `FALSE`, default `FALSE`: assemble a trace across the
  stages. The stages must also have been built with `keep_trace = TRUE`
  to contribute anything.

## Value

An S7 object of class
[Chain](https://statmodels7.github.io/optimizers7/reference/Chain-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).
Its `name` joins the names of the stages with `then`, and its
`criterion`, `maxit` and `max_eval` are copied from the last stage.

## Each stage keeps its own settings

A stage carries its own stopping rule and its own budgets, which is why
a chain is built from optimizers instead of taking arguments: a coarse
rule and a small budget for the exploration, a tight rule for the
descent, written as
`chain(sa(maxit = 200), newton(criterion = crit_grad(1e-12)))`.

## The result

The point, the value, the gradient and `iterations` are those of the
last stage, where the run ended. So is `converged`: a chain has
converged when the method that finished it converged. An earlier stage
that exhausts its budget, the ordinary end of a global search, does not
affect the flag.

The evaluation counts are summed over the stages, so a chain can be
compared with a single optimizer on total work. `message` joins the
messages of the stages, and each message is prefixed with its stage
number.

A trace is assembled when both the chain and the stage were built with
`keep_trace = TRUE`: if the chain was built with `keep_trace = FALSE`,
the result has no trace, and if a stage was built with
`keep_trace = FALSE`, that stage contributes no rows. The trace carries
a `stage` column. Stages need not report the same columns, since a
derivative-free method has a `stationarity` column where a descent has
`gnorm`; when they differ, only the trace of the last stage is kept.

An error raised by a stage propagates, since a method that cannot run on
the objective indicates a problem with the objective. A stage that runs
without converging passes its point on, which the first stage of a chain
usually does. The rule of every stage is checked by
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
before any of them runs, so a chain whose third stage cannot evaluate
its criterion fails before running the first two.

`chain(x)` with a single stage reports the run of that stage: the point,
value, gradient, counts, iterations and verdict are those of the stage.
The `message` carries the prefix `stage 1: `, and `optimizer` is the
chain.

## See also

[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
for the other wrapper,
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) for
the search with which a chain usually begins,
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
for the rule check.

## Examples

``` r
chain(sa(maxit = 10), bfgs())
#> <optimizer> simulated annealing (uniform) then BFGS
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 500, evaluations Inf
#>   settings  : stages = <list>

# Rastrigin has many local minima (121 in the box [-5.12, 5.12]^2), so the
# minimum that a descent method finds depends on its start. The search
# chooses the basin, and the descent finds the minimum inside it.
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

# The counts are the work of the whole chain; the iteration count is that
# of the last stage.
set.seed(3)
res <- minimize(chain(sa(maxit = 20), bfgs()), rastrigin, c(3.5, -2.5))
unlist(res@counts)
#>    f    g    h 
#> 1250    0    0 
res@iterations
#> [1] 4

# A chain of one stage reports the point of that stage.
identical(minimize(chain(bfgs()), rastrigin, c(3.5, -2.5), gr = rg)@par,
          minimize(bfgs(), rastrigin, c(3.5, -2.5), gr = rg)@par)
#> [1] TRUE

# A stage that cannot evaluate its own rule is caught before anything runs.
try(minimize(chain(bfgs(), nelder_mead(criterion = crit_grad())),
             rastrigin, c(1, 1), gr = rg))
#> Error : The stopping rule needs gradient, which nelder-mead does not provide.
#>   Choose a criterion that this optimizer can evaluate, or a method that provides it.
```
