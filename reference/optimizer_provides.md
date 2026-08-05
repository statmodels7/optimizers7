# What an Optimizer Can Offer a Stopping Rule

The names of the `state` components an optimizer is able to fill in, so
that
[`check_criterion`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
can refuse a rule it could never satisfy.

## Usage

``` r
optimizer_provides(optimizer)
```

## Arguments

- optimizer:

  An
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Value

A character vector. The names the shipped criteria read are `"gradient"`
and `"stationarity"`.

## Details

Gradient-based methods provide a gradient. The derivative-free ones
provide a stationarity measure instead, since no single derivative they
could report goes to zero at a solution.

Every optimizer evaluates the objective, so there is no token for that
and rules reading it are never refused. A user-defined method inherits
the default, which claims a gradient; if that is not true of it, say so,
because the refusal machinery relies on this declaration being accurate.

## See also

[`check_criterion`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md),
[`crit_needs`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)

## Examples

``` r
optimizer_provides(bfgs())
#> [1] "gradient"
optimizer_provides(nelder_mead())
#> [1] "stationarity"
```
