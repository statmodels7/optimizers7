# Check a Line Search's Resolution

A single non-negative finite number, zero meaning the question is not
asked, or a function of no arguments returning one.

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

The function form is for an objective whose resolution MOVES. It is
asked once per invocation of the search rather than per trial, so it
costs one call an iteration, and it is what lets a caller whose
objective settles as it goes report the resolution of the current point
instead of the reading from the worst-located point of the run.

## See also

[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
