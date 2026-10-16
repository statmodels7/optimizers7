# What a Derivative-Free Method Can Offer a Stopping Rule

[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
and
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
report `"stationarity"` and no gradient, since they do not compute a
gradient. A rule that reads a gradient is rejected when the run starts;
without this check it would test `NULL` at every iteration and never
fire.

## Arguments

- optimizer:

  A `NelderMead`, `Compass` or `Brent` object.

## Value

The character vector `"stationarity"`.

## Details

The measure differs by method, and
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
reads the measure of the method: for
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
the size of the simplex, the largest distance from the best vertex to
another vertex, so the tolerance is on the parameter scale; for
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
the poll size \\\Delta\\, the quantity in which Torczon's theorem is
stated; and for
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
the half-width of the bracket.

## Examples

``` r
optimizer_provides(nelder_mead())
#> [1] "stationarity"
optimizer_provides(compass())
#> [1] "stationarity"

# A gradient rule is therefore rejected, and the error names the method.
try(minimize(nelder_mead(criterion = crit_grad()),
             function(p) sum(p^2), c(1, 1)))
#> Error : The stopping rule needs gradient, which nelder-mead does not provide.
#>   Choose a criterion that this optimizer can evaluate, or a method that provides it.
```
