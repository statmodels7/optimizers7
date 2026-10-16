# Warn When a Supplied Gradient Does Not Belong to the Objective

Compares the directional derivative of `fn` along the gradient direction
at `par`, obtained from one central difference, with the rate predicted
by the supplied `gr`. A gradient computed from a different model than
the objective, such as one that holds a predictor fixed where the
objective varies it or that reads an outdated copy of the data, makes
the two disagree by a large factor. Without this check the effect is a
line-search failure at the first iteration, with no indication of its
cause.

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

The check costs one call to `gr` and two to `fn`, runs once per call of
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
and is skipped when the comparison cannot be made: for an objective that
is not a function, a gradient that is not finite at `par`, an objective
that is not finite at the probe points, and a gradient too small to be
distinguished from the truncation error of the difference. The last case
arises when a caller starts at the optimum, where a warning would be
wrong.

The tolerance is loose, a relative disagreement above one half, so that
the error of the difference or a subgradient of a non-smooth objective
does not trigger it. It is meant to detect a wrong function, not an
error in the eighth digit. A disagreement at the first step is checked
again at a step ten times smaller and, if it persists, at a step a
hundred times smaller, and the warning is given only if it persists at
both: the step is scaled by the length of `par`, so a large coordinate
can make it too long for a strongly curved direction, and the truncation
error of the difference shrinks with the step, whereas the disagreement
of a wrong gradient does not. Each repetition costs two more calls to
`fn`.

Setting `options(optimizers7.check_gradient = FALSE)` disables it.
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
disables it for all its inner runs, since the generic has already
checked the gradient at the starting value of the caller.
