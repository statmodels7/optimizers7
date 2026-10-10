# What a Criterion Needs From the Iteration

Returns the names of the `state` components that a criterion requires,
so that an algorithm can reject a rule that it cannot evaluate.

## Usage

``` r
crit_needs(criterion)
```

## Arguments

- criterion:

  A
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object.

## Value

A character vector of `state` component names, possibly empty.

## Details

A derivative-free method computes no gradient, so
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
given to such a method would test `NULL` at every iteration and never
stop the run.
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
compares the value of this generic with the value of
[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
and, when the run starts, signals an error that names the missing
component and the method.

The base method returns
[`character()`](https://rdrr.io/r/base/character.html), so a rule that
reads only the objective needs no method here: every optimizer evaluates
the objective, and such a rule is accepted by every method. The shipped
rules declare two names, `"gradient"` and `"stationarity"`.

## See also

[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
for the rule itself,
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
for the rejection this feeds,
[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
for the other half of the comparison.

## Examples

``` r
# A rule reading only the objective declares nothing; a gradient rule
# declares the one component it reads.
crit_needs(crit_abs_obj())
#> character(0)
crit_needs(crit_grad())
#> [1] "gradient"

# A combination declares the union, so a method that lacks any one of
# them rejects it.
crit_needs(crit_any(crit_grad(), crit_stationary()))
#> [1] "gradient"     "stationarity"
```
