# Validate a Positive Whole Number

Checks that the value is a single positive whole number no larger than
`.Machine$integer.max`. It is used for `max_step`, which counts trial
points, and for the `memory` of
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md),
which counts secant pairs. The count is passed to the compiled code as
an integer, so a fractional value is rejected (`memory = 2.5` would
otherwise run at a memory of 2), and so is an infinite one, which
[`as.integer()`](https://rdrr.io/r/base/integer.html) would turn into
`NA`.

## Usage

``` r
check_count(v, nm)
```

## Arguments

- v:

  The value.

- nm:

  Its name, for the message.

## Value

Invisibly `TRUE`. Raises an error naming `nm` otherwise.

## See also

[`check_whole()`](https://statmodels7.github.io/optimizers7/reference/check_whole.md)
for a count that may be zero.
