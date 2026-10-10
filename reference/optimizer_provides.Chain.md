# What a Chain Can Offer a Stopping Rule

Returns the answer of the last stage. The stopping rule of the chain is
that of the last stage, so that stage has to be able to evaluate it, and
a chain that begins with a derivative-free search can still carry a
gradient rule.

## Arguments

- optimizer:

  A `Chain` object.

## Value

A character vector of `state` component names, those of the last stage.

## Examples

``` r
# A search then a descent: the rule reads the gradient of the descent.
optimizer_provides(chain(sa(), bfgs()))
#> [1] "gradient"

# The other order ends with nelder_mead(), which offers the stationarity
# measure instead.
optimizer_provides(chain(bfgs(), nelder_mead()))
#> [1] "stationarity"
```
