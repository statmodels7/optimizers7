# S7 Class for a Combination of Criteria

The class
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
and
[`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)
both build, and the two methods it implements. Its
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
evaluates every rule it holds and reduces with
[`any()`](https://rdrr.io/r/base/any.html) or
[`all()`](https://rdrr.io/r/base/all.html) according to `how`; its
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)
is the union of what they need, so a combination containing a gradient
rule is refused by a derivative-free method exactly as the bare rule
would be.

## Usage

``` r
CritCombine(label = character(0), criteria = list(), how = character(0))
```

## Arguments

- criteria:

  A list of
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  objects.

- how:

  Either `"any"` or `"all"`.

## Value

An S7 object of class `CritCombine` inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
carrying `label`, `criteria` and `how`.

## Details

A `CritCombine` is itself a
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
so combinations nest and the label nests with them:
`crit_any(crit_all(a, b), c)` reads `a and b or c`. Every rule is
evaluated at every call, [`any()`](https://rdrr.io/r/base/any.html) and
[`all()`](https://rdrr.io/r/base/all.html) taking the whole vector
rather than short-circuiting, which costs nothing worth counting against
an objective evaluation.

## See also

[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
and
[`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)
for the constructors,
[`combine_criteria()`](https://statmodels7.github.io/optimizers7/reference/combine_criteria.md)
for the shared body.

## Examples

``` r
# The needs are the union, so this is refused by a simplex method.
crit_needs(crit_any(crit_grad(), crit_stationary()))
#> [1] "gradient"     "stationarity"

st <- list(iter = 3, f_new = 1, f_old = 2, x_new = 1, x_old = 1,
           gradient = c(1e-9, -2e-9))
c(any = crit_met(crit_any(crit_grad(1e-8), crit_never()), st),
  all = crit_met(crit_all(crit_grad(1e-8), crit_never()), st))
#>   any   all 
#>  TRUE FALSE 

# Combinations nest, and so does the label.
crit_any(crit_all(crit_grad(), crit_abs_par()), crit_never())@label
#> [1] "gradient (max-norm) < 1e-06 and |dx| < 1e-08 or iteration budget"
```
