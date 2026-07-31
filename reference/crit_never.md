# Never Stop Early

The rule that never fires, so a run ends only when it exhausts its
iteration budget.

## Usage

``` r
crit_never()
```

## Value

A
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
object.

## Details

This is not a placeholder. For a stochastic method there is often
nothing left to test: every quantity a convergence rule could look at —
the objective, the gradient — is a noisy estimate drawn from whichever
observations happened to be sampled, and a tolerance applied to one of
those measures the noise rather than the progress. Such a run is meant
to be governed by its budget, and saying so with an object is better
than leaving a real criterion in place that quietly never fires.

A run that ends this way reports `converged = FALSE`, which is the
truth: the budget ran out, and nothing checked whether the answer was
any good. It is the same discipline everywhere else in the package —
convergence is what a rule confirmed, never what the run merely stopped
doing.

## See also

[`adam`](https://statmodels7.github.io/optimizers7/reference/adam.md),
[`crit_grad`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)

## Examples

``` r
crit_never()
#> <criterion> iteration budget
```
