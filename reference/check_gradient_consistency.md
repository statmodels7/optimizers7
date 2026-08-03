# Warn When a Supplied Gradient Does Not Belong to the Objective

Compares the directional derivative of `fn` along the gradient direction
at `par`, obtained from one central difference, with the rate the
supplied `gr` predicts. A gradient computed from a different model than
the objective – a fixed predictor where the objective uses the
parameter, a stale copy of the data – makes the two disagree grossly,
and the symptom downstream is otherwise a mute line-search failure at
the first iteration.

## Usage

``` r
check_gradient_consistency(fn, gr, par)
```

## Arguments

- fn:

  The objective.

- gr:

  The supplied gradient.

- par:

  The starting point, numeric.

## Value

Invisibly `TRUE`; warns on a gross mismatch.

## Details

The check costs one call to `gr` and two to `fn`, runs once per
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
call, and is skipped whenever it cannot be decisive: a non-function
objective, a zero or non-finite gradient at `par`, an objective that is
not finite at the probe points. The tolerance is deliberately loose – a
relative disagreement above one half – so that finite-difference error
or a subgradient of a non-smooth objective does not trip it; it exists
to catch the wrong function, not the eighth digit.

Setting `options(optimizers7.check_gradient = FALSE)` disables it, and
[`multistart`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
disables it for its inner runs after the first.
