# S7 Class for the Result of an Optimization

What
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
and
[`maximize()`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
return: the point reached and the value there, the evaluation counts,
whether a stopping rule was satisfied and which one, and enough of the
run to diagnose it when none was. The optimizer itself is kept on the
object, so a run can be repeated from what it returned.

## Usage

``` r
optimizer_result(
  par = integer(0),
  value = integer(0),
  gradient = NULL,
  counts = NULL,
  iterations = integer(0),
  converged = logical(0),
  criterion_met = character(0),
  message = character(0),
  trace = NULL,
  optimizer = NULL,
  elapsed = integer(0),
  seed = NULL
)
```

## Arguments

- par:

  The minimizer reached, a numeric vector of the same length as the
  starting point.

- value:

  The objective at `par`, a single number. For a run started by
  [`maximize()`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
  this is the value of the original objective, with the sign already
  undone.

- gradient:

  The gradient at `par`, a numeric vector, or `NULL` for a method that
  computes none.
  [`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
  reports the proximal gradient mapping here.

- counts:

  A named integer vector of length three, `f`, `g` and `h`: evaluations
  of the objective, the gradient and the Hessian. A gradient obtained by
  differencing is counted in `f`, one per coordinate per side.

- iterations:

  The number of iterations performed, a single integer.

- converged:

  A single logical; see Details.

- criterion_met:

  A single string: the label of the rule that fired, or
  `iteration budget reached` or `evaluation budget exhausted` when the
  run ended without one.

- message:

  A single string, empty when there is nothing to report. Carries notes
  such as `gradient obtained by finite differences`.

- trace:

  The iteration path as a data frame, or `NULL` when the optimizer was
  built with `keep_trace = FALSE`, which is the default.

- optimizer:

  The
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  that produced the result, kept whole so the run can be repeated or
  restarted from where it stopped.

- elapsed:

  Wall-clock seconds, a single number. Measured to the platform's clock
  resolution, so a fast run can read exactly `0`.

- seed:

  The state of the random number generator when the run began, an
  integer vector, or `NULL` for a deterministic method.

## Value

An S7 object of class `optimizer_result` carrying the twelve properties
above. Built by the methods of
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md);
a caller receives one rather than constructing it.

## `converged` means a rule fired

`converged` is `TRUE` only when the stopping rule was satisfied. Running
out of iterations or of evaluations leaves it `FALSE` and puts
`iteration budget reached` or `evaluation budget exhausted` into
`criterion_met`. Reporting a budget as a success is the commonest defect
in a hand-written optimization loop; it turns a failure into a wrong
answer that nothing downstream can detect.

## The trace, and the column that changes name

With `keep_trace = TRUE` on the optimizer, `trace` is a data frame with
one row per iteration and five columns: `iteration`, `value`, the
quantity the stopping rule is watching, `step` and `safeguard`. The
third column is `gnorm` for a method that computes a gradient and
`stationarity` for one that does not, so code reading a trace should ask
for the column by position or check `names(trace)`.

`safeguard` names the repair the algorithm applied to its own step, or
`"none"`. The names differ by method and are worth reading:
`step shortened` and `step adjusted` for the line-search methods,
`cg restart` when conjugacy is lost, `bb curvature reset` when a secant
pair reports none, and `reflect`, `expand`, `contract in` and
`contract out` for the simplex.
[`summary()`](https://rdrr.io/r/base/summary.html) tabulates them.

## Repeating a stochastic run

`seed` holds `.Random.seed` as it stood when the run began, and is
filled in only by the methods that draw:
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md), a
resampling
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md),
a
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
generating its own starts. Assigning it back reproduces the run exactly:

    a <- minimize(sa(maxit = 300), f, c(-1.2, 1))
    assign(".Random.seed", a@seed, globalenv())
    b <- minimize(sa(maxit = 300), f, c(-1.2, 1))
    identical(a@par, b@par)   # TRUE

For a deterministic method `seed` is `NULL`.

## See also

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
[`print.optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/print.optimizer_result.md)
and
[`summary.optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/summary.optimizer_result.md)
for the two views of it,
[`plot.optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/plot.optimizer_result.md)
for the trace.

## Examples

``` r
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))

res <- minimize(bfgs(keep_trace = TRUE), f, c(-1.2, 1), gr = g)
res@par
#> [1] 0.9999995 0.9999989
res@converged
#> [1] TRUE
res@criterion_met
#> [1] "gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08"
res@counts
#>  f  g  h 
#> 49 39  0 

# A budget stops the run and leaves converged FALSE, with the reason.
short <- minimize(bfgs(maxit = 3), f, c(-1.2, 1), gr = g)
c(short@converged, short@criterion_met)
#> [1] "FALSE"                    "iteration budget reached"

# The trace names its third column after what the rule watches.
names(res@trace)
#> [1] "iteration" "value"     "gnorm"     "step"      "safeguard"
names(minimize(nelder_mead(keep_trace = TRUE), f, c(-1.2, 1))@trace)
#> [1] "iteration"    "value"        "stationarity" "step"         "safeguard"   

# The optimizer travels with the result, so the run can be repeated.
again <- minimize(res@optimizer, f, c(-1.2, 1), gr = g)
all.equal(again@par, res@par)
#> [1] TRUE
```
