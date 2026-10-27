# S7 Class for the Ridge Modification

This class holds the smallest shift tried by the ridge modification.
[`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md)
builds it.

## Usage

``` r
ModifyRidge(label = character(0), floor = integer(0))
```

## Arguments

- label:

  The label of every
  [`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md).

- floor:

  The smallest shift, as in
  [`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md).

## Value

An S7 object of class `ModifyRidge` inheriting from
[`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md),
with the property `floor` beside `label`.

## See also

[`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md)
for the constructor.
