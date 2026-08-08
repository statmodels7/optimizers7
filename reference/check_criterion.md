# Reject a Stopping Rule an Optimizer Cannot Evaluate

Compares what the optimizer's criterion reads against what the optimizer
can supply, and stops if the rule asks for something absent.

## Usage

``` r
check_criterion(optimizer)
```

## Arguments

- optimizer:

  An
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Value

Invisibly `TRUE`; raises an error otherwise.

## Details

Call this at the top of a user-defined
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
method. A rule needing a gradient, handed to a method that computes
none, would sit testing `NULL` at every iteration and never fire; the
run would then end on its iteration budget and report a reason nowhere
near the truth. Rejecting it here, by name, is what
[`check_optimizer`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
tests for.

What an optimizer can supply is declared by
[`optimizer_provides`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md),
whose default claims a gradient. Override that when the method computes
none.

## See also

[`optimizer_provides`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md),
[`crit_needs`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)

## Examples

``` r
check_criterion(bfgs())
try(check_criterion(nelder_mead(criterion = crit_grad())))
#> Error : The stopping rule needs gradient, which nelder-mead does not provide.
#>   Choose a criterion this optimizer can evaluate, or a method that provides it.
```
