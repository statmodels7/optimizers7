# Validate a Tolerance

Checks that the value is a single positive number, so that every
criterion constructor rejects the same nonsense in the same words. Zero
is refused as well as a negative: a rule with a tolerance of zero can
never fire.

## Usage

``` r
check_tol(tol)
```

## Arguments

- tol:

  The value supplied.

## Value

Invisibly `TRUE`. Raises an error naming `tol` otherwise.

## Details

It is also called from outside this file, by
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
for its proximity weights and by
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) for
its step bounds, so the message names `tol` rather than the caller's own
argument. A reader who wrote `bb(alpha0 = 0)` sees a message about
`tol`.
