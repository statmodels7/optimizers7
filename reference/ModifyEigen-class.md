# S7 Class for the Eigenvalue Modification

This class holds the lower bound of the eigenvalues of the modified
Hessian.
[`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md)
builds it.

## Usage

``` r
ModifyEigen(label = character(0), floor = integer(0))
```

## Arguments

- label:

  The label of every
  [`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md).

- floor:

  The lower bound, as in
  [`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md).

## Value

An S7 object of class `ModifyEigen` inheriting from
[`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md),
with the property `floor` beside `label`.

## See also

[`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md)
for the constructor.
