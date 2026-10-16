# What the Bundle Method Can Offer a Stopping Rule

Reports `"stationarity"` and not `"gradient"`, although the method
evaluates subgradients and reports their aggregate. A rule built on
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
would test the subgradient at the iterate, which does not tend to zero:
near the minimum of \\\lvert x \rvert\\ every evaluated subgradient has
norm 1, so the rule would not fire when the run reached the minimizer.

## Arguments

- optimizer:

  A `Bundle` object.

## Value

The character vector `"stationarity"`.

## Details

[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
reads the optimality estimate \\\lVert p \rVert^{2} + \alpha\\ instead,
which tends to zero, and is the default rule of this method.
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
also reads this declaration before its second check, which therefore
does not require the `gradient` reported by a bundle run to equal the
gradient of the objective.

## Examples

``` r
optimizer_provides(bundle())
#> [1] "stationarity"

# A gradient rule is therefore rejected when the run starts.
f <- function(p) sum(abs(p - c(1, -2)))
try(minimize(bundle(criterion = crit_grad()), f, c(0, 0),
             gr = function(p) sign(p - c(1, -2))))
#> Error : The stopping rule needs gradient, which proximal bundle does not provide.
#>   Choose a criterion that this optimizer can evaluate, or a method that provides it.
```
