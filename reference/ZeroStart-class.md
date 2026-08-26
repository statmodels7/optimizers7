# S7 Class for the Zero Starter

A starter that produces a vector of zeros on the unconstrained scale.
Built by
[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md).
It carries `npar` alone, adding no property of its own to the abstract
[starter](https://statmodels7.github.io/optimizers7/reference/starter-class.md)
class.

## Usage

``` r
ZeroStart(npar = NULL)
```

## Arguments

- npar:

  The number of parameters, an integer, or `NULL`.

## Value

An S7 object of class `ZeroStart` inheriting from
[starter](https://statmodels7.github.io/optimizers7/reference/starter-class.md).

## See also

[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
for the constructor,
[`starting_values.ZeroStart()`](https://statmodels7.github.io/optimizers7/reference/starting_values.ZeroStart.md)
for what it produces.
