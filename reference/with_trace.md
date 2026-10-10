# Rebuild an Optimizer With the Trace Switched On

Returns a copy of the optimizer with `keep_trace = TRUE`, keeping its
class and every other setting.
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
uses it for check 6, which tests that the trace is well formed and
therefore needs a run that records one.

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
The method for
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
also switches on the trace of every stage, because the trace of a chain
is assembled from the traces of its stages. A
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
needs no method of its own: its trace is the table of starts, which the
wrapper builds from the results of the inner runs.
