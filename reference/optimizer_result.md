# S7 Class for the Result of an Optimisation

What
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
returns: the answer, how it was reached, and enough of the run to
diagnose it when it was not reached.

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

  The minimiser.

- value:

  The objective there.

- gradient:

  The gradient there, or `NULL` if the method computes none.

- counts:

  Evaluations of the objective and of the gradient.

- iterations:

  Iterations performed.

- converged:

  Logical; see Details.

- criterion_met:

  Which rule ended the run.

- message:

  A human-readable account.

- trace:

  The iteration path, or `NULL`.

- optimizer:

  The
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  that produced this, kept so the run can be repeated exactly.

- elapsed:

  Seconds.

- seed:

  The state of the random number generator when the run began, for a
  method that draws any, and `NULL` otherwise.

## Value

An S7 object of class `optimizer_result`.

## Details

`converged` is `TRUE` only when the stopping rule was satisfied. It is
**never** `TRUE` because the iteration budget ran out — that is the
commonest defect in hand-written optimisation loops, and it turns a
failure into a silently wrong answer.

`trace`, present when the optimiser was built with `keep_trace = TRUE`,
records for each iteration the objective, the quantity the criterion is
watching, the step actually taken, and the name of any safeguard that
fired. That last column is what turns "it did not converge" into a
diagnosis.

`seed` is filled in by the methods that draw random numbers – a `"mads"`
poll, a subsampling
[`adam`](https://statmodels7.github.io/optimizers7/reference/adam.md), a
[`multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
generating its own starts. Assigning it back with
`assign(".Random.seed", res@seed, globalenv())` reproduces the run
exactly. A stochastic method that cannot be repeated is very hard to
debug, and remembering to call
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand is not
something anyone does until the second time they need it.

## See also

[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

## Examples

``` r
res <- minimize(bfgs(), function(p) sum((p - c(1, 2))^2), c(0, 0),
                gr = function(p) 2 * (p - c(1, 2)))
res@par
#> [1] 1 2
res@converged
#> [1] TRUE
res@criterion_met
#> [1] "gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)"
res@counts
#> f g h 
#> 3 2 0 
```
