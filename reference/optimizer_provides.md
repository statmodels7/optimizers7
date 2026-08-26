# What an Optimizer Can Offer a Stopping Rule

The names of the `state` components an optimizer is able to fill in, so
that
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
can reject a rule it could never satisfy.

## Usage

``` r
optimizer_provides(optimizer)
```

## Arguments

- optimizer:

  An
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Value

A character vector of `state` component names. The two the shipped
criteria read are `"gradient"` and `"stationarity"`.

## Details

The gradient-based methods provide `"gradient"`. The derivative-free
ones provide `"stationarity"` instead, no single derivative they could
report going to zero at a solution:
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
offers the simplex diameter,
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
the poll size,
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md)
Corana's termination measure and
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
its optimality estimate.

Every optimizer evaluates the objective, so there is no token for that
and a rule reading the objective is never rejected.
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
answers `"gradient"` and means the proximal gradient mapping, which
vanishes at a stationary point of the whole objective.

The **default method** answers `"gradient"`, so an optimizer written
outside the package claims one unless it says otherwise. Say otherwise
if it is untrue: the rejection machinery relies on this declaration
being accurate, and a method that claims a gradient it does not compute
will accept a rule that can never fire.

## See also

[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
for the rejection this feeds,
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)
for the other half of the comparison,
[`optimizer_bounded()`](https://statmodels7.github.io/optimizers7/reference/optimizer_bounded.md)
for the package's other declaration generic.

## Examples

``` r
# The two answers, across the shipped methods.
vapply(list(bfgs(), newton(), cg(), bb(), gd(), adam()),
       optimizer_provides, "")
#> [1] "gradient" "gradient" "gradient" "gradient" "gradient" "gradient"
vapply(list(nelder_mead(), compass(), sa(), bundle()),
       optimizer_provides, "")
#> [1] "stationarity" "stationarity" "stationarity" "stationarity"

# A wrapper answers for whichever run reports the result.
optimizer_provides(chain(sa(), bfgs()))
#> [1] "gradient"
optimizer_provides(multistart(nelder_mead()))
#> [1] "stationarity"
```
