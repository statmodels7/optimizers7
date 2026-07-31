# S7 Class for a Combination of Criteria

The class
[`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
and
[`crit_all`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)
instantiate.

## Usage

``` r
CritCombine(label = character(0), criteria = list(), how = character(0))
```

## Arguments

- criteria:

  A list of
  [`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  objects.

- how:

  Either `"any"` or `"all"`.

## Value

An S7 object inheriting from
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

## See also

[`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md),
[`crit_all`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)
