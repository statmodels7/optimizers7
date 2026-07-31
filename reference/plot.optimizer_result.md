# Plot Method for an Optimisation Result

The objective against iteration, with any iteration at which a safeguard
fired marked. Requires `keep_trace = TRUE`.

## Arguments

- x:

  An
  [`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

- ...:

  Passed to [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

No return value; called for the plot.

## Examples

``` r
res <- minimize(gd(keep_trace = TRUE),
                function(p) sum((p - 1:2)^2), c(0, 0))
plot(res)

```
