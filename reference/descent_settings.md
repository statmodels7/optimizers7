# The Settings of a Descent Method for the Compiled Loop

`descent_settings()` returns the arguments that a descent method passes
to the compiled loop for a problem with `npar` parameters, without
running it, for another package that runs the loop through the headers
in `inst/include/optimizers7` (`descent_loop()` in `descent_loop.h`).

## Usage

``` r
descent_settings(optimizer, npar)

descent_capture
```

## Format

`descent_capture`: an environment whose flag `on` makes
[`run_descent()`](https://statmodels7.github.io/optimizers7/reference/run_descent.md)
return its settings instead of running.

## Arguments

- optimizer:

  An
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

- npar:

  The number of parameters.

## Value

A list, or `NULL`.

## Details

The settings are produced by the method's own
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
method, which stops before the run, so they are the ones a call of
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
with the same optimizer would use: the description of the direction
(`method`), the line search as
[`line_search_spec()`](https://statmodels7.github.io/optimizers7/reference/line_search_spec.md)
describes it (`line_search`), `maxit`, `max_eval`
(`.Machine$integer.max` for an unlimited budget), `verbose`, `refresh`,
`keep_trace`, `step`, and the stopping rule as
[`criterion_spec()`](https://statmodels7.github.io/optimizers7/reference/criterion_spec.md)
describes it (`rule_spec`, `NULL` where the rule has no compiled form).
The result is `NULL` for an optimizer whose method does not run the
descent loop, and for
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md),
whose direction reads a third derivative supplied as an R function.

## See also

[`criterion_spec()`](https://statmodels7.github.io/optimizers7/reference/criterion_spec.md),
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

## Examples

``` r
str(descent_settings(lbfgs(), 2)$method)
#> List of 4
#>  $ type          : chr "lbfgs"
#>  $ memory        : int 10
#>  $ curvature_tol : num 1e-10
#>  $ stop_decrement: num 0.001
descent_settings(nelder_mead(), 2)
#> NULL
```
