# Run Optimizers One After Another

Wraps several optimizers into one, each starting where the previous
finished.

## Usage

``` r
chain(..., verbose = FALSE, keep_trace = FALSE)
```

## Arguments

- ...:

  Two or more optimizers, in the order they should run. A single one is
  accepted.

- verbose:

  Report which stage is running? Defaults to `FALSE`.

- keep_trace:

  Store the path of every stage? Defaults to `FALSE`.

## Value

An S7 object of class `Chain`, inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Details

The composition a global search needs: `chain(sa(), lbfgs())` explores
first and then descends from wherever the exploration left off, and
neither method has to know about the other. It is the second wrapper of
this shape after
[`multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.md),
and the two compose – `multistart(chain(sa(), lbfgs()))` is a legal
optimizer.

Each stage carries its OWN criterion and its own budgets, which is the
point of a chain rather than an argument: a coarse rule and a small
budget for the exploration, a tight rule for the descent.

**What the result reports.** The point and the value are the LAST
stage's, since that is where the run ended, and so is `converged`: a
chain has converged when the method that finished it says so, and an
earlier stage exhausting its own budget is the ordinary way a global
search ends rather than a failure of the whole. Evaluations and
iterations are summed over the stages, and the trace, when kept, carries
a `stage` column.

A stage that raises propagates, since a method that cannot run on the
objective is a fact about the objective. A stage that runs without
converging passes its point on, which is what the first stage of a chain
usually does.

`chain(x)` with a single stage is that stage's run, reported through the
chain: the wrapper is not a special case to be avoided.

## See also

[`multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.md),
[`sa`](https://statmodels7.github.io/optimizers7/reference/sa.md)

## Examples

``` r
chain(sa(maxit = 10), bfgs())
#> <optimizer> simulated annealing (uniform) then BFGS
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 500, evaluations Inf
#>   settings  : stages = <list>

# a multimodal objective: the search says which basin, the descent finishes
# the job inside it
rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
set.seed(3)
minimize(sa(maxit = 20), rastrigin, c(3.5, -2.5))@value
#> [1] 1.479952
set.seed(3)
minimize(chain(sa(maxit = 20), bfgs()), rastrigin, c(3.5, -2.5))@value
#> [1] 0.9949591
```
