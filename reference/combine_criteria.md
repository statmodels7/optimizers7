# The Shared Body of the Two Combinators

Validates the arguments and builds the combined criterion, so that
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
and
[`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)
reject the same nonsense in the same words: an empty call and an
argument that is not a criterion.

## Usage

``` r
combine_criteria(dots, how)
```

## Arguments

- dots:

  A list of
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  objects.

- how:

  Either `"any"` or `"all"`.

## Value

An S7 object of class
[CritCombine](https://statmodels7.github.io/optimizers7/reference/CritCombine-class.md),
inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

## Details

The label is the sub-labels joined by `or` or `and` according to `how`,
which is why a nested combination reads as one sentence.
