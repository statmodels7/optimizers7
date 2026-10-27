# S7 Class for Hessian Modifications

This class is the abstract parent of the objects that repair an
indefinite or nearly singular Hessian before a Newton step.
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
hold one in their `modify` argument. The two shipped subclasses are
built by
[`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md)
and
[`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md);
a new modification is a subclass of this class with a method for
[`hessian_modify()`](https://statmodels7.github.io/optimizers7/reference/hessian_modify.md).

## Usage

``` r
hessian_modification(label = character(0))
```

## Arguments

- label:

  A short character label, shown when the optimizer that contains the
  modification is printed.

## Value

An S7 object of class `hessian_modification`, carrying `label`. The
class is abstract, so every value is an object of a subclass.

## Details

A modification is used only at an iterate where the Cholesky
factorization of the Hessian fails. Where the factorization succeeds,
the Newton step is computed from the Hessian itself and no modification
is involved. See
[`hessian_modify()`](https://statmodels7.github.io/optimizers7/reference/hessian_modify.md)
for the requirements on a subclass and for an example of a user-defined
modification, and
[modifications](https://statmodels7.github.io/optimizers7/reference/modifications.md)
for a comparison of the shipped ones.

## See also

[`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md),
[`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md),
[`hessian_modify()`](https://statmodels7.github.io/optimizers7/reference/hessian_modify.md),
[modifications](https://statmodels7.github.io/optimizers7/reference/modifications.md),
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)

## Examples

``` r
# The class is abstract: use a constructor or a subclass.
try(hessian_modification(label = "mine"))
#> Error in new_object(S7_object(), label = label) : 
#>   Can't construct an object from abstract class <hessian_modification>

modify_eigen()
#> <hessian_modification> eigenvalue floor (floor = 1e-08)
modify_ridge(floor = 1e-6)
#> <hessian_modification> ridge (floor = 1e-06)
```
