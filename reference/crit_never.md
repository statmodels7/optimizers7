# Never Stop Early

Builds the rule that never fires, so a run ends only when it exhausts a
budget and reports `converged = FALSE`. This is
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)'s
default and is the honest arrangement for a stochastic method, where
every quantity a convergence rule could read is an estimate.

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

This is not a placeholder. For a stochastic method there is often
nothing left to test: every quantity a convergence rule could look at,
the objective and the gradient alike, is a noisy estimate drawn from
whichever observations happened to be sampled, and a tolerance applied
to one of those measures the noise. Such a run is meant to be governed
by its budget, and saying so with an object beats leaving a real
criterion in place that quietly never fires.

A run that ends this way reports `converged = FALSE`, which is the
truth: the budget ran out, and nothing checked whether the answer was
any good. That is the package's rule everywhere. Convergence is what a
stopping rule confirmed, never what the run merely stopped doing.

## See also

[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
for the method that defaults to it,
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
for the rule to pass instead when the objective is exact.

## Examples

``` r
crit_never()
#> <criterion> iteration budget

# It fires at no state at all, however good.
crit_met(crit_never(), list(f_new = 0, f_old = 0, x_new = 1, x_old = 1,
                            gradient = c(0, 0)))
#> [1] FALSE

# A run carrying it ends on its budget and says so.
r <- minimize(adam(maxit = 50), function(p) sum(p^2), c(1, 1),
              gr = function(p) 2 * p)
c(r@converged, r@criterion_met)
#> [1] "FALSE"                    "iteration budget reached"
```
