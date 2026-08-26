# Validate a Constant in the Unit Interval

Checks that the value is a single number **strictly** inside \\(0, 1)\\.
Both endpoints are refused: `c1 = 0` asks for no decrease at all and
`c1 = 1` asks for the whole decrease the linear model predicts, which a
curved objective cannot supply.

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
