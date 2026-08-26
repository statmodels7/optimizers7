# S7 Class for the Uniform Starter

A starter that draws each coordinate independently from a uniform on the
unconstrained scale. Built by
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md).
It adds `min` and `max` to the `npar` the abstract
[starter](https://statmodels7.github.io/optimizers7/reference/starter-class.md)
class carries; either may be one number or one per parameter.

## Usage

``` r
UniformStart(npar = NULL, min = integer(0), max = integer(0))
```

## Arguments

- npar:

  The number of parameters, an integer, or `NULL`.

- min, max:

  The range drawn from, in unconstrained units.

## Value

An S7 object of class `UniformStart` inheriting from
[starter](https://statmodels7.github.io/optimizers7/reference/starter-class.md).

## See also

[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md)
for the constructor,
[`starting_values.UniformStart()`](https://statmodels7.github.io/optimizers7/reference/starting_values.UniformStart.md)
for the draw.
