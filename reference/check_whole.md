# Validate a Non-Negative Whole Number

Checks that the value is a single whole number, zero included, no larger
than `.Machine$integer.max`. Used where a count says how many of
something to tolerate, so that zero means that none is tolerated: the
`memory` of
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md),
where zero makes the reference the current value and the search
identical to
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md),
and the `max_skip` of
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
where zero resets the approximation on the first skipped update.

A negative value is rejected because the compiled code treats it as
zero: it compares `max_skip` with `>=` and tests `memory > 0`. An
infinite value is rejected because
[`as.integer()`](https://rdrr.io/r/base/integer.html) would turn it into
`NA`.

## Usage

``` r
check_whole(v, nm)
```

## Arguments

- v:

  The value.

- nm:

  Its name, for the message.

## Value

Invisibly `TRUE`. Raises an error naming `nm` otherwise.

## See also

[`check_count()`](https://statmodels7.github.io/optimizers7/reference/check_count.md)
for a count that has to be at least one.
