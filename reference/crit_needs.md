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
  [`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object.

## Value

A character vector, possibly empty.

## Details

A derivative-free method has no gradient, so
[`crit_grad`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
would sit there testing `NULL` at every iteration and quietly never stop
the run. Refusing it at construction is the same discipline that makes
`check_link()` in linkfunctions7 report a numerical derivative order as
numerical rather than as passed: a check that cannot be evaluated must
say so.

## Examples

``` r
crit_needs(crit_grad(1e-8))
#> [1] "gradient"
crit_needs(crit_rel_obj(1e-10))
#> character(0)
```
