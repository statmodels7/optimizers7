# Stop Only When Every Rule Fires

Combines criteria conjunctively, for a run that should not stop until
several independent things agree.

## Usage

``` r
crit_all(...)
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

[`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)

## Examples

``` r
crit_all(crit_grad(1e-6), crit_abs_par(1e-10))
#> <criterion> gradient (max-norm) < 1e-06 and |dx| < 1e-10
```
