# Validate a Decrement Threshold

Checks that `stop_decrement` is a single positive number, or `Inf`, so
that
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md),
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md),
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
reject invalid values with the same message.

## Usage

``` r
check_stop_decrement(x)
```

## Arguments

- x:

  The value supplied.

## Value

Invisibly `TRUE`. Raises an error naming `stop_decrement` otherwise.
