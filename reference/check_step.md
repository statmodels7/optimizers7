# Validate an Initial Step Length

Checks that `step` is a single positive number, so that all six methods
taking one reject the same nonsense in the same words.

## Usage

``` r
check_step(step)
```

## Arguments

- step:

  The value supplied.

## Value

Invisibly `TRUE`. Raises an error naming `step` otherwise.
