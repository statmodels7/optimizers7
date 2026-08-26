# Plot Method for an Optimization Result

Draws the objective against the iteration number as a line, with a
filled point at every iteration where a safeguard fired and a legend
saying so. The method's name is the title. Reading the marks against the
curve shows where the algorithm was in trouble and whether the objective
moved when it was.

## Arguments

- x:

  An
  [`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).
  `type`, `lwd`, `las`, `xlab`, `ylab` and `main` are already set and
  passing them again is an error, as it is for any duplicated argument.

## Value

`NULL`, invisibly. Called for the plot.

## Details

The plot is built from the `trace`, so the optimizer must have been
created with `keep_trace = TRUE`. Without one the method stops with a
message naming that argument, since an empty plot would say nothing.

The vertical axis is the objective on its own scale; a run whose
objective falls over several orders of magnitude is better read with
`log = "y"`, which passes through to
[`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html) like
any other argument.

## Examples

``` r
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))

# The marked iterations are the ones where a step had to be repaired.
plot(minimize(newton(keep_trace = TRUE), f, c(-1.2, 1), gr = g), log = "y")


# Without a trace there is nothing to draw, and the method says which
# argument was missing.
try(plot(minimize(newton(), f, c(-1.2, 1), gr = g)))
#> Error : No trace to plot; build the optimizer with keep_trace = TRUE.
```
