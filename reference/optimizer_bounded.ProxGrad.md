# The Proximal Gradient Method Does Not Take Box Bounds

Returns `FALSE`, the only shipped method that does. A box constraint is
itself a non-smooth term, expressed by the projection onto the box, and
this method already has a slot for such a term: `prox`. Offering bounds
beside it would be a second and conflicting route to the same thing, and
the two operators would have to be composed by somebody.

## Arguments

- optimizer:

  A `ProxGrad` object.

## Value

`FALSE`.

## Details

[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
consults this before running its bounds check, so a method that answers
`FALSE` is not failed for refusing a box it never promised. The
composition a caller writes instead is
`prox(v, t) = pmin(pmax(prox_penalty(v, t), lower), upper)`, which
imposes both terms exactly.

## Examples

``` r
pg <- prox_grad(prox = function(v, t) v, g = function(b) 0)
optimizer_bounded(pg)
#> [1] FALSE
optimizer_bounded(bfgs())
#> [1] TRUE
```
