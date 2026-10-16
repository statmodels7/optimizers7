# S7 Class for Optimization Algorithms

This class is the abstract parent of every algorithm in the package. An
optimizer object carries the algorithm and every setting that the
algorithm uses: its stopping rule, its two budgets, what it reports
during the run and what it keeps. It is built once and can be reused,
printed, stored beside a result and passed around, so that all the
settings of a run are stored in the object and not in the call that
started it.

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
  into the result. Set by each constructor; a caller does not supply it.

- criterion:

  The stopping rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object. Anything else raises an error naming
  [`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
  as an example. An algorithm also rejects a rule that it cannot
  evaluate, through
  [`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md).

- maxit:

  Maximum iterations, a single number at least 1 and finite. `Inf` is
  rejected, because a run without a last resort can run forever on an
  objective that never satisfies its rule.

- max_eval:

  Maximum evaluations of the objective, a single number at least 1.
  `Inf` is the default and is allowed, so the budget is off unless a
  finite value is given. A run stopped here reports `converged = FALSE`,
  because the budget ended it and not the rule.

- verbose:

  `TRUE` or `FALSE`; whether to report progress during the run. `NA` is
  rejected.

- refresh:

  Report every `refresh` iterations when `verbose` is `TRUE`. It is a
  single non-negative number; `0` is allowed and reports only the final
  summary.

- keep_trace:

  `TRUE` or `FALSE`; whether to store the iteration path in the `trace`
  field of the result. It is off by default, since the path costs one
  row per iteration.

## Value

An S7 object of class `optimizer`, with properties `name` (character),
`criterion` (a
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)),
`maxit`, `max_eval` and `refresh` (numeric) and `verbose` and
`keep_trace` (logical). The class is abstract, so it is never returned
directly: every value is an object of one of its subclasses.

## The seven properties of every algorithm

Each algorithm is a subclass that adds its own settings: a Hessian
repair for
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
a curvature tolerance for
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
the two moment decay rates for
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md).
The seven properties below are common to all of them, so a caller can
replace one optimizer with another and leave the rest of the call
unchanged. [`print()`](https://rdrr.io/r/base/print.html) shows the
stopping rule and the two budgets, and then the settings of the
algorithm on a `settings` line; `verbose`, `refresh` and `keep_trace`
are not shown.

## The iteration and evaluation budgets

`maxit` bounds the iterations and `max_eval` bounds the evaluations of
the objective, and the two differ whenever a line search is expensive:
one iteration may spend many evaluations, and the number of backtracks
does not count towards `maxit`. A run stopped by `max_eval` reports
`converged = FALSE`, with `criterion_met` equal to
`evaluation budget reached` (see the examples).

`maxit` must be finite, because it is the stop of last resort;
`max_eval` may be `Inf`, which is its default.

## See also

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
for the run,
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
for the stopping rule,
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
for the checks that a new subclass has to pass.

## Examples

``` r
# The class is abstract, so it cannot be built: use one of the constructors.
try(optimizer(name = "mine"))
#> Error in new_object(S7_object(), name = name, criterion = criterion, maxit = maxit,  : 
#>   Can't construct an object from abstract class <optimizer>

# Every algorithm has the same seven properties, and each one has properties
# of its own.
shared <- c("name", "criterion", "maxit", "max_eval", "verbose",
            "refresh", "keep_trace")
stopifnot(all(shared %in% names(S7::props(nelder_mead()))))
setdiff(names(S7::props(newton())), shared)
#> [1] "step"        "line_search" "hessian_mod" "floor"       "max_length" 
#> [6] "typical"    
setdiff(names(S7::props(adam())), shared)
#> [1] "alpha"   "beta1"   "beta2"   "eps"     "decay"   "amsgrad"

# Printing shows the rule and the budgets before the settings of the
# algorithm.
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
