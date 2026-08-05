# Rebuild an Optimizer With a Different Stopping Rule

Replaces the criterion, and for a wrapper replaces the one that will
actually be consulted.

## Usage

``` r
with_criterion(optimizer, criterion)
```

## Arguments

- optimizer:

  The
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

- criterion:

  The new rule.

## Value

An optimizer of the same class.

## Details

The distinction matters.
[`multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
carries a criterion only so that printing it tells the truth; the rule
that is evaluated belongs to the optimizer inside. Setting the outer one
and expecting a different run is the sort of thing that makes a check
pass while testing nothing.
