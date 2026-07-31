# What the Bundle Method Can Offer a Stopping Rule

The objective and the predicted decrease, but not a gradient.

## Arguments

- optimizer:

  A `Bundle` object.

## Value

A character vector.

## Details

It evaluates subgradients and reports their aggregate, but it does not
offer `"gradient"`, and the omission is deliberate.
[`crit_grad`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
would then test a quantity that does not go to zero — at the minimum of
\\\lvert x \rvert\\ every subgradient has norm 1 — so the rule would
never fire while sitting on the answer.
[`crit_stationary`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
tests the predicted decrease instead, which does go to zero.
