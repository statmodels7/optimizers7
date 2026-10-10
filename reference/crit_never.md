# The Rule That Never Fires

Builds the rule that never fires, so a run ends when it exhausts a
budget and reports `converged = FALSE`. It is the default of
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md),
a method for stochastic objectives, on which every quantity that a
stopping rule could read is an estimate.

## Usage

``` r
crit_never()
```

## Value

An S7 object of class
[CritNever](https://statmodels7.github.io/optimizers7/reference/CritNever-class.md),
inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

## Details

On a stochastic objective the value and the gradient are noisy estimates
computed from the observations that were sampled, and a tolerance
applied to one of them tests the noise. Such a run is governed by its
budget, and this rule states that explicitly.

A run that ends on its budget reports `converged = FALSE`, because no
rule tested the answer. A descent method can also stop before its
budget, when its line search accepts no step, and it then reports
`converged = FALSE` as well. In the package `converged` is `TRUE` only
when a stopping rule fired, when a line search with a declared
resolution established that no step decreases the objective by more than
that resolution, or when
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
reduced its bracket to the resolution of the iterate.

## See also

[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
for the method that defaults to it,
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
for the rule to pass instead when the objective is exact.

## Examples

``` r
crit_never()
#> <criterion> iteration budget

# It fires at no state, even at an exact minimum.
crit_met(crit_never(), list(f_new = 0, f_old = 0, x_new = 1, x_old = 1,
                            gradient = c(0, 0)))
#> [1] FALSE

# A run with this rule ends on its budget, and the result reports it.
r <- minimize(adam(maxit = 50), function(p) sum(p^2), c(1, 1),
              gr = function(p) 2 * p)
c(r@converged, r@criterion_met)
#> [1] "FALSE"                    "iteration budget reached"
```
