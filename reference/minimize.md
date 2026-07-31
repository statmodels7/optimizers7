# Minimise a Function

The entry point of the package. Everything here minimises; see
[`maximize`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
for the other direction.

## Usage

``` r
minimize(
  optimizer,
  fn,
  par,
  gr = NULL,
  he = NULL,
  lower = -Inf,
  upper = Inf,
  ...
)
```

## Arguments

- optimizer:

  An
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  object, carrying the algorithm and its settings.

- fn:

  The objective: a function of the parameter vector returning a single
  number.

- par:

  A numeric vector of starting values.

- gr:

  An optional gradient function. Ignored when `fn` carries its own
  gradient.

- he:

  An optional Hessian function. Methods that do not use one ignore it,
  so calling code need not branch on the algorithm.

- lower, upper:

  Box constraints. Numeric, of length one — applying to every parameter
  — or one value per parameter. Default `-Inf` and `Inf`, which is no
  constraint at all. See Details.

- ...:

  Passed to methods.

## Value

An
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

## Details

The generic dispatches on `optimizer` alone, so each algorithm is
written once. The objective is normalised separately by
[`as_objective`](https://statmodels7.github.io/optimizers7/reference/as_objective.md),
which dispatches on `fn`, so a caller with its own kind of objective
registers one method there and every algorithm accepts it.

A gradient that is not supplied is computed by central finite
differences. The result records which derivatives were supplied and
which were differenced, so a run is never silently less exact than it
appears.

`lower` and `upper` are two vectors rather than a list of pairs, which
is what [`optim`](https://rdrr.io/r/stats/optim.html) and
[`nlminb`](https://rdrr.io/r/stats/nlminb.html) take, and what lets
`lower = 0` say "every parameter is positive" without writing out one
pair per coefficient. Recycling is length one or one per parameter and
nothing between, since anything else is far more likely to be a mistake
than a request.

**Bounds are removed, not enforced.** Each bounded coordinate is
reparametrised onto the whole real line — a shifted log for one-sided
bounds, a scaled logit for two — and the optimiser runs unconstrained in
the new variable. Every point it proposes is admissible by construction,
so there is no rejection step and no boundary for a line search to trip
over, and any method works with bounds without knowing about them.

The limitation is worth knowing before it surprises you: **an optimum
lying on a bound cannot be reached**. Getting there requires the
transformed variable to run to infinity, so the optimiser marches off,
improves by less and less, and stops on a budget at a point merely close
to the bound. For the statistical use this exists to serve — a positive
variance, a probability inside the unit interval — the optimum is
interior and this never arises. For a genuine box-constrained problem
with active constraints at the solution, an active-set method is the
right tool and this is not one.

## See also

[`maximize`](https://statmodels7.github.io/optimizers7/reference/maximize.md),
[`gd`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)

## Examples

``` r
# a quadratic, with the gradient supplied
minimize(gd(), function(p) sum((p - c(1, 2))^2),
         par = c(0, 0), gr = function(p) 2 * (p - c(1, 2)))
#> <optimizer_result> gradient descent
#>   value      : 0
#>   par        : 1 2
#>   iterations : 1   evaluations: f 3, g 2
#>   converged  : yes (gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative))

# and without, so the gradient is differenced
minimize(gd(), function(p) sum((p - c(1, 2))^2), c(0, 0))
#> <optimizer_result> gradient descent
#>   value      : 2.80957e-22
#>   par        : 1 2
#>   iterations : 1   evaluations: f 11, g 0
#>   converged  : yes (gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative))
#>   note       : gradient obtained by finite differences

# one bound for every parameter: a scale that must stay positive
minimize(bfgs(), function(p) sum((p - c(1, 2))^2), c(0.5, 0.5), lower = 0)
#> <optimizer_result> BFGS
#>   value      : 4.68179e-20
#>   par        : 1 2
#>   iterations : 7   evaluations: f 42, g 0
#>   converged  : yes (gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative))
#>   note       : gradient obtained by finite differences

# or one per parameter. The unconstrained minimum is at (1, 2), so the
# second coordinate is pushed against its ceiling of 1.
minimize(bfgs(), function(p) sum((p - c(1, 2))^2), c(0.5, 0.5),
         lower = c(0, 0), upper = c(5, 1))
#> <optimizer_result> BFGS
#>   value      : 1
#>   par        : 1 1
#>   iterations : 30   evaluations: f 161, g 0
#>   converged  : yes (gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative))
#>   note       : gradient obtained by finite differences
```
