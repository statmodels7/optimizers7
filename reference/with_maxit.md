# Rebuild an Optimizer With a Different Iteration Budget

Returns a copy of the optimizer with `maxit` replaced, keeping its class
and every other setting.
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
uses it for check 3, which gives an optimizer a single iteration and
tests that it does not report convergence, and which receives an
optimizer of any class.

## Usage

``` r
with_maxit(optimizer, maxit)
```

## Arguments

- optimizer:

  The
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  to copy.

- maxit:

  The new budget, a single finite number at least 1.

## Value

An optimizer of the same class as `optimizer`.

## Details

The default method sets the property with
[`S7::set_props()`](https://rconsortium.github.io/S7/reference/props.html),
which is correct for an optimizer whose own budget governs the run. A
wrapper needs a method of its own so that the change reaches the
optimizer inside;
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
and
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
have one.
