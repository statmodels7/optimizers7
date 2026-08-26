# What a Chain Can Offer a Stopping Rule

Reports whatever the last stage offers. The chain's stopping rule is the
last stage's, so it is that stage which has to be able to evaluate it,
and a chain opening with a derivative-free search can still carry a
gradient rule.

## Arguments

- optimizer:

  A `Chain` object.

## Value

A character vector of `state` component names, the last stage's.

## Examples

``` r
# A search then a descent: the descent's gradient is what the rule reads.
optimizer_provides(chain(sa(), bfgs()))
#> [1] "gradient"

# The other order offers only the simplex-free measure.
optimizer_provides(chain(bfgs(), nelder_mead()))
#> [1] "stationarity"
```
