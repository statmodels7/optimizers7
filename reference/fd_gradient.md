# A Central-Difference Gradient

Computes the finite-difference gradient that the proximal loop uses when
the caller supplies none. The step is scaled to each coordinate.

## Usage

``` r
fd_gradient(fn, x)
```

## Arguments

- fn:

  A function of the parameter vector.

- x:

  The point.

## Value

A numeric vector.
