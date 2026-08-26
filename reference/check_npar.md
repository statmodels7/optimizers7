# Validate a Declared Parameter Count

Checks that `npar` is a single positive whole number, or `NULL`, and
returns it as an integer. Called by both starter constructors, so
`start_zeros(0)` and `start_runif(npar = 2.5)` are refused in the same
words at the point they are written.

## Usage

``` r
check_npar(npar)
```

## Arguments

- npar:

  `NULL` or a positive whole number.

## Value

`NULL` when `npar` is `NULL`, otherwise the value as an integer. Raises
an error otherwise.
