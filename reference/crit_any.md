# Stop When Any of Several Rules Fires

Combines criteria disjunctively. This is the usual arrangement: a run
should end as soon as any reasonable rule is satisfied.

## Usage

``` r
crit_any(...)
```

## Arguments

- ...:

  [`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  objects.

## Value

A
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
object, so combinations nest.

## See also

[`crit_all`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)

## Examples

``` r
crit_any(crit_grad(1e-8), crit_rel_obj(1e-12))
#> <criterion> gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
```
