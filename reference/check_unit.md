# Validate a Constant in the Unit Interval

Checks that the value is a single number **strictly** inside \\(0, 1)\\.
Both endpoints are rejected: `c1 = 0` requires no decrease at all, and
`c1 = 1` requires the whole decrease that the linear model predicts,
which a strictly convex objective cannot supply.

## Usage

``` r
check_unit(v, nm)
```

## Arguments

- v:

  The value.

- nm:

  Its name, for the message.

## Value

Invisibly `TRUE`. Raises an error naming `nm` otherwise.
