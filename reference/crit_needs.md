# What a Criterion Needs From the Iteration

The names of the `state` components a criterion requires, so that an
algorithm can reject a rule it cannot evaluate instead of accepting one
that never fires.

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

A derivative-free method has no gradient, so
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
handed to one would sit there testing `NULL` at every iteration and
quietly never stop the run.
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
compares what this reports against what
[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
reports and refuses the pair when the run starts, naming both.

The base method returns
[`character()`](https://rdrr.io/r/base/character.html), so a rule
reading only the objective needs no method here: every optimizer
evaluates the objective and no rule reading it is ever refused. The two
names the shipped rules declare are `"gradient"` and `"stationarity"`.

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

# A combination declares the union, so a method missing any one of them
# is refused.
crit_needs(crit_any(crit_grad(), crit_stationary()))
#> [1] "gradient"     "stationarity"
```
