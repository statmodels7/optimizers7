# What a Derivative-Free Method Can Offer a Stopping Rule

Both derivative-free methods report `"stationarity"` and no gradient,
having none to report. A rule reading a gradient is refused when the run
starts, so it cannot sit there testing `NULL` at every iteration and
never firing.

## Arguments

- optimizer:

  A `NelderMead` or `Compass` object.

## Value

The character vector `"stationarity"`.

## Details

The measure differs by method and
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
reads whichever is offered: for
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
the **diameter of the simplex**, so the tolerance is on the parameter
scale, and for
[`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md)
the **poll size** \\\Delta\\, the quantity Torczon's theorem is stated
in.

## Examples

``` r
optimizer_provides(nelder_mead())
#> [1] "stationarity"
optimizer_provides(compass())
#> [1] "stationarity"

# So a gradient rule is refused, naming the method.
try(minimize(nelder_mead(criterion = crit_grad()),
             function(p) sum(p^2), c(1, 1)))
#> Error : The stopping rule needs gradient, which nelder-mead does not provide.
#>   Choose a criterion this optimizer can evaluate, or a method that provides it.
```
