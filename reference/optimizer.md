# S7 Class for Optimization Algorithms

The abstract parent of every algorithm in the package. An optimizer
object carries the algorithm and every setting the algorithm obeys: its
stopping rule, its two budgets, what it reports as it goes and what it
keeps. It is built once and can be reused, printed, stored beside a
result and passed around, so nothing about a run is hidden inside the
call that started it.

## Usage

``` r
optimizer(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0)
)
```

## Arguments

- name:

  A short character name, shown when the object is printed and carried
  into the result. Set by each constructor; a caller has no reason to
  supply it.

- criterion:

  The stopping rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object. Anything else raises an error naming
  [`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
  as an example. An algorithm also rejects a rule it cannot evaluate,
  through
  [`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md).

- maxit:

  Maximum iterations, a single number at least 1 and **finite**. `Inf`
  is refused: a run with no last resort can hang on an objective that
  never satisfies its rule.

- max_eval:

  Maximum evaluations of the objective, a single number at least 1.
  `Inf` is the default and is allowed, so the budget is off unless it is
  asked for. A run stopped here reports `converged = FALSE`, the budget
  having ended it rather than the rule.

- verbose:

  `TRUE` or `FALSE`; whether to report progress as the run goes. `NA` is
  refused.

- refresh:

  Report every `refresh` iterations when `verbose` is `TRUE`. A single
  non-negative number; `0` is allowed and reports only the final
  summary.

- keep_trace:

  `TRUE` or `FALSE`; whether to store the iteration path in the result's
  `trace` field. Off by default, the path costing one row per iteration.

## Value

An S7 object of class `optimizer`, with properties `name` (character),
`criterion` (a
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)),
`maxit`, `max_eval` and `refresh` (numeric) and `verbose` and
`keep_trace` (logical). The class is abstract, so it is never returned
directly: every value is an object of one of its subclasses.

## The seven properties every algorithm has

Each algorithm is a subclass adding the settings that are its own:
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
a Hessian repair,
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
a curvature tolerance,
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
its three decay rates. The seven below are common to all of them, so a
caller can replace one optimizer with another and leave every other line
alone. [`print()`](https://rdrr.io/r/base/print.html) shows the shared
settings first and the algorithm's own on a `settings` line.

## The two budgets are different budgets

`maxit` bounds progress and `max_eval` bounds work, and they diverge
whenever a line search is expensive: one iteration may spend many
evaluations, and a search that backtracks thirty times is invisible to
`maxit`. On Rosenbrock from the customary start,
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
with no cap takes 35 iterations and 210 evaluations; capped at
`max_eval = 40` it stops after 7 iterations and 42 evaluations,
reporting `converged = FALSE` and `evaluation budget exhausted`.

`maxit` must be finite, because it is the stop of last resort;
`max_eval` may be `Inf`, which is its default.

## See also

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
for the run,
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
for the stopping rule,
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
for what a subclass of your own has to promise.

## Examples

``` r
# Abstract, so it cannot be built: use one of the constructors.
try(optimizer(name = "mine"))
#> Error in new_object(S7_object(), name = name, criterion = criterion, maxit = maxit,  : 
#>   Can't construct an object from abstract class <optimizer>

# The same seven properties on every algorithm, plus each one's own.
shared <- c("name", "criterion", "maxit", "max_eval", "verbose",
            "refresh", "keep_trace")
stopifnot(all(shared %in% names(S7::props(nelder_mead()))))
setdiff(names(S7::props(newton())), shared)
#> [1] "step"        "line_search" "hessian_mod" "floor"      
setdiff(names(S7::props(adam())), shared)
#> [1] "alpha"   "beta1"   "beta2"   "eps"     "decay"   "amsgrad"

# Printing shows the rule and the budgets before anything algorithm-specific.
bfgs()
#> <optimizer> BFGS
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 500, evaluations Inf
#>   settings  : step = 1, line_search = strong Wolfe (c1 = 1e-04, c2 = 0.9), curv_tol = 1e-10, max_skip = 5

# The two budgets bound different things.
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
a <- minimize(bfgs(maxit = 1000), f, c(-1.2, 1))
b <- minimize(bfgs(maxit = 1000, max_eval = 40), f, c(-1.2, 1))
c(iterations = a@iterations, evaluations = a@counts[["f"]])
#>  iterations evaluations 
#>          35         210 
c(iterations = b@iterations, evaluations = b@counts[["f"]])
#>  iterations evaluations 
#>           7          42 
b@message
#> [1] "evaluation budget exhausted; gradient obtained by finite differences"
```
