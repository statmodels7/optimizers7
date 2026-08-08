# A Central-Difference Gradient

The finite-difference gradient the R-level proximal loop uses when the
caller supplies none, with the step scaled to each coordinate.

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
