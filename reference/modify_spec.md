# Describe a Modification to the Compiled Loop

Returns the part of the `method` list of
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
that describes the Hessian modification: `modify_type` (`"eigen"`,
`"ridge"` or `"custom"`), `floor`, and `modify_fn`, `NULL` for the two
shipped modifications and otherwise a function of the Hessian that calls
[`hessian_modify()`](https://statmodels7.github.io/optimizers7/reference/hessian_modify.md).

## Usage

``` r
modify_spec(modify)
```

## Arguments

- modify:

  A
  [`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md)
  object.

## Value

A list of three.
