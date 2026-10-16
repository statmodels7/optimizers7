# Whether a Chain Takes Box Bounds

Returns `TRUE` when every stage takes bounds and `FALSE` otherwise. The
bounds are passed to all of them, so a single stage that does not take
them makes the whole chain unable to take them.

## Arguments

- optimizer:

  A `Chain` object.

## Value

A single logical.

## Details

[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
is the only shipped method that returns `FALSE`: it takes its constraint
inside the proximal operator, where the constraint combines with the
penalty. A chain that contains it therefore does not take bounds,
whichever other stages it has.

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
