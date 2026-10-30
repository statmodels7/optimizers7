# Plot Method for an Optimization Result

[`autoplot()`](https://statmodels7.github.io/optimizers7/reference/autoplot.md)
returns a graph of the objective against the iteration number, drawn as
a line, with a point at every iteration where a safeguard fired;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints that
graph. The title is the name of the method. The points show the
iterations at which the algorithm repaired its step and how the
objective moved there. A legend below the panel names the line and the
points.

## Usage

``` r
# S3 method for class 'optimizer_result'
plot(x, legend = TRUE, ...)

autoplot_optimizer_result(object, legend = TRUE, ...)
```

## Arguments

- x, object:

  An
  [`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

- legend:

  Whether the graph draws the legend below the panel, `TRUE` (the
  default) or `FALSE`.

- ...:

  Unused; an argument passed here signals an error.

## Value

[`autoplot()`](https://statmodels7.github.io/optimizers7/reference/autoplot.md)
returns a ggplot2 object;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) prints the same
object and returns it invisibly.

## Details

The graph is built from the `trace`, so the optimizer must have been
created with `keep_trace = TRUE`. Without a trace the method stops with
an error that names that argument.

The vertical axis is the objective on its own scale. For a run whose
objective falls over several orders of magnitude, a logarithmic axis is
more readable, and it is added to the returned graph with
`+ ggplot2::scale_y_log10()`.

The graph is drawn with ggplot2, which is suggested; the methods signal
an error that names it when it is not installed.

## Examples

``` r
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))
res <- minimize(newton(keep_trace = TRUE), f, c(-1.2, 1), gr = g)

# The points mark the iterations at which a step was repaired.
p <- plot(res)


# The objective falls by twenty orders of magnitude, so a logarithmic
# axis shows the last iterations.
p + ggplot2::scale_y_log10()


# Without a trace there is nothing to draw, and the error names the
# missing argument.
try(plot(minimize(newton(), f, c(-1.2, 1), gr = g)))
#> Error : Plotting needs a trace; build the optimizer with keep_trace = TRUE.
```
