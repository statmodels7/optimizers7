# Whether a Multistart Accepts Box Bounds

Returns the answer of the inner optimizer, since every start is a run of
that optimizer with the same bounds. A multistart over
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
therefore does not accept box bounds.

## Arguments

- optimizer:

  A `MultiStart` object.

## Value

A single logical.
