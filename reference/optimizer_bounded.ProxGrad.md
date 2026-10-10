# The Proximal Gradient Method Does Not Take Box Bounds

Returns `FALSE`; it is the only shipped method that does. A box
constraint is itself a non-smooth term, expressed by the projection onto
the box, and this method already receives such a term through `prox`.
Bounds beside it would be a second, conflicting way of imposing a
constraint, and the two operators would have to be composed.

## Arguments

- optimizer:

  A `ProxGrad` object.

## Value

`FALSE`.

## Details

[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
reads this declaration before its bounds check, so a method that returns
`FALSE` is not failed for rejecting bounds. For a penalty that is convex
and separable coordinate by coordinate, such as the lasso or the ridge,
the composition
`prox(v, t) = pmin(pmax(prox_penalty(v, t), lower), upper)` is the
proximal operator of the penalty plus the box. For a nonconvex or a
non-separable penalty the composition is not the operator of the sum.

## Examples

``` r
pg <- prox_grad(prox = function(v, t) v, g = function(b) 0)
optimizer_bounded(pg)
#> [1] FALSE
optimizer_bounded(bfgs())
#> [1] TRUE
```
