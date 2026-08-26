# Is This a Starter?

`TRUE` when `x` inherits from the abstract
[starter](https://statmodels7.github.io/optimizers7/reference/starter-class.md)
class, `FALSE` otherwise.
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
asks it to decide whether `par` is a vector to be used as given or an
object to be resolved into one, and the question has a name so that the
test is written once.

## Usage

``` r
is_starter(x)
```

## Arguments

- x:

  Any object.

## Value

A single logical.
