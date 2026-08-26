# Print Method for Criteria

Shows a criterion in one line, as `<criterion>` and its label. The label
is what the result reports in `criterion_met`, so printing a rule shows
exactly the string a converged run will carry.

## Arguments

- x:

  A
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object.

- ...:

  Unused.

## Value

`x`, invisibly. Called for the output.

## Examples

``` r
crit_grad()
#> <criterion> gradient (max-norm) < 1e-06
crit_stationary(1e-10)
#> <criterion> stationarity < 1e-10

# A combination prints as one sentence, and nests.
crit_any(crit_grad(), crit_rel_obj())
#> <criterion> gradient (max-norm) < 1e-06 or |df| < 1e-12 (relative)
crit_any(crit_all(crit_grad(), crit_abs_par()), crit_never())
#> <criterion> gradient (max-norm) < 1e-06 and |dx| < 1e-08 or iteration budget

# And it is the string the result reports.
r <- minimize(bfgs(criterion = crit_grad()), function(p) sum(p^2), c(1, 1),
              gr = function(p) 2 * p)
r@criterion_met
#> [1] "gradient (max-norm) < 1e-06"
```
