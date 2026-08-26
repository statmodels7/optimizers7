# Recycle a Length-One Vector, and Reject Any Other Mismatch

Returns `v` at length `n`: a single value is repeated, a value already
of length `n` is passed through as a double, and anything else raises an
error naming the argument and both lengths. R's own recycling is
deliberately not used, since it is silent whenever the shorter length
divides the longer, and a partial range is far likelier to be a mistake
than a request.

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

  The argument's name, for the message.

## Value

A numeric vector of length `n`.
