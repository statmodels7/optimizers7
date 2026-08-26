# Rebuild an Optimizer With the Trace Switched On

Returns a copy of the optimizer with `keep_trace = TRUE`, keeping its
class and every other setting.
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
needs it for check 6, which asks whether the trace is well formed, and a
trace has to be asked for.

## Usage

``` r
with_trace(optimizer)
```

## Arguments

- optimizer:

  The
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  to copy.

## Value

An optimizer of the same class as `optimizer`, with `keep_trace` `TRUE`.

## Details

The default method sets the property with
[`S7::set_props()`](https://rconsortium.github.io/S7/reference/props.html).
A wrapper needs a method of its own so that the inner optimizer records
a path too;
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
has one.
