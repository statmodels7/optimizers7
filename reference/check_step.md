# Validate an Initial Step Length

Checks that `step` is a single positive number, so that every
constructor that takes one rejects invalid values with the same message.

## Usage

``` r
check_step(step)
```

## Arguments

- step:

  The value supplied.

## Value

Invisibly `TRUE`. Raises an error naming `step` otherwise.
