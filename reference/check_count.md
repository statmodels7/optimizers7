# Validate a Positive Whole Number

Checks that the value is a single positive integer. Used for `max_step`,
which counts trials, and for
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)'s
`memory`, which counts secant pairs. A fractional value is a mistake in
both places and is refused: the kernel receives the count through
[`as.integer()`](https://rdrr.io/r/base/integer.html), so `memory = 2.5`
would otherwise run silently at a memory of 2.

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
for a count whose zero is meaningful.
