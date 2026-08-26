# Rebuild an Optimizer With a Different Stopping Rule

Returns a copy of the optimizer with its criterion replaced, keeping its
class and every other setting. For a wrapper it replaces the rule that
is actually consulted, which is not always the one on the outside.

## Usage

``` r
with_criterion(optimizer, criterion)
```

## Arguments

- optimizer:

  The
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  to copy.

- criterion:

  The new rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object.

## Value

An optimizer of the same class as `optimizer`.

## Details

[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
carries a criterion so that printing it tells the truth; the rule the
run evaluates belongs to the optimizer inside. Setting the outer one
alone changes the printing and nothing else, which is exactly the sort
of thing that makes a check pass while testing nothing:

    ms <- multistart(bfgs(), n = 3)
    with_criterion(ms, crit_abs_obj(1e-4))@optimizer@criterion@label
    # "|df| < 1e-04"  (the inner rule changed too)

    S7::set_props(ms, criterion = crit_abs_obj(1e-4))@optimizer@criterion@label
    # "gradient (max-norm) < 1e-06 or ..."  (unchanged)

[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
has a method of its own for the same reason, its reported rule being the
last stage's.
