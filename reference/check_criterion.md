# Reject a Stopping Rule That an Optimizer Cannot Evaluate

Compares the quantities that the criterion of the optimizer reads with
the quantities that the optimizer supplies, and signals an error if the
rule needs a quantity that the optimizer does not supply.

## Usage

``` r
check_criterion(optimizer)
```

## Arguments

- optimizer:

  An
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Value

Invisibly `TRUE`; raises an error otherwise.

## Details

Call this at the top of a user-defined
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
method. A rule that needs a gradient, used with a method that computes
none, would test `NULL` at every iteration and never fire; the run would
then end on its iteration budget and report a failure unrelated to the
state of the run.
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
tests that a method rejects such a rule.

The quantities that an optimizer supplies are declared by
[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md),
whose default method returns `"gradient"`. Override it when the method
computes no gradient.

## See also

[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md),
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)

## Examples

``` r
check_criterion(bfgs())
try(check_criterion(nelder_mead(criterion = crit_grad())))
#> Error : The stopping rule needs gradient, which nelder-mead does not provide.
#>   Choose a criterion that this optimizer can evaluate, or a method that provides it.
```
