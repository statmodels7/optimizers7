# Validate a Non-Negative Whole Number

Checks that the value is a single whole number, zero included. Used
where a count says how many of something to tolerate, so that zero says
*tolerate none* and is a legitimate setting:
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)'s
`memory`, where zero makes the reference the current value and the
search is
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
exactly, and
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)'s
`max_skip`, where zero resets the approximation on the first skipped
update.

A negative value is refused because the kernel compares against it with
`>=`, so any negative setting behaves as zero without saying so.

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
