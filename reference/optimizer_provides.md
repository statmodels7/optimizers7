# What an Optimizer Can Offer a Stopping Rule

Returns the names of the `state` components that an optimizer supplies,
so that
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
can reject a rule that the optimizer cannot evaluate.

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

The gradient-based methods provide `"gradient"`.
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md),
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md),
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) and
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
provide `"stationarity"` instead:
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
offers the largest distance from the best vertex to another vertex of
the simplex,
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
the poll size,
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
the half-width of its bracket,
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md)
Corana's termination measure and
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
its optimality estimate.
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
evaluates subgradients but does not provide `"gradient"`, since a
subgradient does not vanish at a kink.

Every optimizer evaluates the objective, so there is no name for it, and
a rule that reads the objective is accepted by every method.
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
returns `"gradient"`, and its gradient is the proximal gradient mapping,
which vanishes at a stationary point of the whole objective.

The **default method** returns `"gradient"`, so an optimizer written
outside the package declares a gradient unless it has a method of its
own. The rejection relies on this declaration being accurate: a method
that declares a gradient that it does not compute accepts a rule that
can never fire.

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

# A wrapper returns the answer of the run that reports the result.
optimizer_provides(chain(sa(), bfgs()))
#> [1] "gradient"
optimizer_provides(multistart(nelder_mead()))
#> [1] "stationarity"
```
