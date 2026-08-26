# What the Bundle Method Can Offer a Stopping Rule

Reports `"stationarity"` and withholds `"gradient"`, though the method
evaluates subgradients and reports their aggregate. The omission is
deliberate:
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
would test a quantity that never goes to zero, since at the minimum of
\\\lvert x \rvert\\ every subgradient has norm 1, so the rule would sit
on the answer without firing.

## Arguments

- optimizer:

  A `Bundle` object.

## Value

The character vector `"stationarity"`.

## Details

[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
reads the optimality estimate \\\lVert p \rVert^{2} + \alpha\\ instead,
which does go to zero, and is this method's default rule.
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
also consults this function before its second check, so a bundle run is
not asked to prove that what it reports as `gradient` is one.

## Examples

``` r
optimizer_provides(bundle())
#> [1] "stationarity"

# So a gradient rule is refused when the run starts.
f <- function(p) sum(abs(p - c(1, -2)))
try(minimize(bundle(criterion = crit_grad()), f, c(0, 0),
             gr = function(p) sign(p - c(1, -2))))
#> Error : The stopping rule needs gradient, which proximal bundle does not provide.
#>   Choose a criterion this optimizer can evaluate, or a method that provides it.
```
