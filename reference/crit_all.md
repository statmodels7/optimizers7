# Stop Only When Every Rule Fires

Combines criteria conjunctively, for a run that should stop only when
several rules hold at once. A conjunction can only become stronger as
terms are added, so a run that ends under it would have ended under any
one of its terms alone at the same iteration or earlier.

## Usage

``` r
crit_all(...)
```

## Arguments

- ...:

  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  objects. At least one is required; anything that is not a criterion
  raises an error.

## Value

An S7 object of class
[CritCombine](https://statmodels7.github.io/optimizers7/reference/CritCombine-class.md),
inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
so combinations nest.

## Details

It is useful when a single rule is known to fire early:
`crit_all(crit_grad(), crit_abs_par())` requires a stationary point at
which the iterate has also stopped moving, a condition that a run
circling a flat optimum does not satisfy.

The needs of a conjunction are the union of the needs of its rules, so
an optimizer rejects a conjunction that contains a rule that it cannot
evaluate, as it rejects that rule alone.

## See also

[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
for the disjunction,
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
for what a rule is.

## Examples

``` r
crit_all(crit_grad(1e-6), crit_abs_par(1e-10))
#> <criterion> gradient (max-norm) < 1e-06 and |dx| < 1e-10

# Both must hold. Here the gradient rule fires and the parameter rule does
# not, so the conjunction does not.
st <- list(f_new = 1, f_old = 2, x_new = c(1, 2), x_old = c(1, 3),
           gradient = c(1e-9, -2e-9))
c(grad = crit_met(crit_grad(1e-8), st),
  par  = crit_met(crit_abs_par(1e-10), st),
  all  = crit_met(crit_all(crit_grad(1e-8), crit_abs_par(1e-10)), st))
#>  grad   par   all 
#>  TRUE FALSE FALSE 
```
