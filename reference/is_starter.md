# Is This a Starter?

Returns `TRUE` when `x` inherits from the abstract
[starter](https://statmodels7.github.io/optimizers7/reference/starter-class.md)
class, and `FALSE` otherwise.
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
calls it to test whether `par` is a vector to be used as given or an
object to be resolved into one, so that the test is written in one
place.

## Usage

``` r
is_starter(x)
```

## Arguments

- x:

  Any object.

## Value

A single logical.
