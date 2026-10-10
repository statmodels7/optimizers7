# Recycle a Length-One Vector, and Reject Any Other Mismatch

Returns `v` at length `n`: a single value is repeated, a value already
of length `n` is passed through as a double, and any other length raises
an error that names the argument and both lengths. The recycling of R is
not used, since it is silent whenever the shorter length divides the
longer, and a partial range is much more likely to be a mistake than a
request.

## Usage

``` r
recycle_to(v, n, nm)
```

## Arguments

- v:

  A numeric vector.

- n:

  The length wanted.

- nm:

  The name of the argument, for the message.

## Value

A numeric vector of length `n`.
