# Validate the Modification of a Newton Method

Checks that `x` inherits from
[`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md).
A character string signals an error that names the constructor to use.

## Usage

``` r
check_modify(x)
```

## Arguments

- x:

  The value supplied.

## Value

Invisibly `TRUE`. Raises an error naming `modify` otherwise.
