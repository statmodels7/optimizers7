# Check a Line Search's Resolution

Checks that `resolution` is a single non-negative finite number, where
zero disables the test, or a function of no arguments that returns one.

## Usage

``` r
check_resolution(x)
```

## Arguments

- x:

  What the constructor was given.

## Value

`NULL`, invisibly; called for the error.

## Details

The function form is for an objective whose resolution changes during
the run. It is called once per invocation of the search and not once per
trial, so it costs one call per iteration. A caller whose objective
becomes more precise during the run can then return the resolution at
the current point instead of the value measured at the start.

## See also

[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
