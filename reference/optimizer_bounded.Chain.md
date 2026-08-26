# Whether a Chain Takes Box Bounds

`TRUE` when every stage takes bounds, `FALSE` otherwise. The bounds are
passed to all of them, so a single stage that would ignore them leaves
the chain unable to promise the box.

## Arguments

- optimizer:

  A `Chain` object.

## Value

A single logical.

## Details

[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
is the one shipped method that answers `FALSE`: it takes its constraint
inside the proximal operator, where the constraint composes with the
penalty already there. A chain containing it is therefore unbounded
whatever else is in it.

## Examples

``` r
optimizer_bounded(chain(bfgs(), newton()))
#> [1] TRUE

pg <- prox_grad(prox = function(v, t) v, g = function(b) 0)
optimizer_bounded(pg)
#> [1] FALSE
optimizer_bounded(chain(bfgs(), pg))
#> [1] FALSE
```
