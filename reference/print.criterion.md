# Print Method for Criteria

Print Method for Criteria

## Arguments

- x:

  A
  [`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object.

- ...:

  Unused.

## Value

`x`, invisibly.

## Examples

``` r
print(crit_any(crit_grad(), crit_rel_obj()))
#> <criterion> gradient (max-norm) < 1e-06 or |df| < 1e-12 (relative)
```
