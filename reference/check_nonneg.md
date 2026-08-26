# Validate a Non-Negative Threshold

Checks that the value is a single finite number that is not negative.
Used for the curvature threshold `curv_tol` of
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md),
which is compared as \\s'y \> \mathrm{curv\\tol} \\ \lVert s \rVert
\lVert y \rVert\\. Zero is admitted and is the textbook curvature
condition \\s'y \> 0\\; a negative value is refused because it makes the
comparison hold for every pair, which turns the skip protection off
without saying so.

## Usage

``` r
check_nonneg(v, nm)
```

## Arguments

- v:

  The value.

- nm:

  Its name, for the message.

## Value

Invisibly `TRUE`. Raises an error naming `nm` otherwise.
