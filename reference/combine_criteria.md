# The Shared Body of the Two Combinators

Validates the arguments and builds the combined criterion, so that
[`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
and
[`crit_all`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)
refuse the same nonsense in the same words.

## Usage

``` r
combine_criteria(dots, how)
```

## Arguments

- dots:

  A list of
  [`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  objects.

- how:

  Either `"any"` or `"all"`.

## Value

A
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
object.
