# The Forward Bound Transform

Maps a parameter from inside its bounds to the unconstrained scale: the
inverse of
[`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)'s
`h`, and what a starting value goes through before an unconstrained run.

## Usage

``` r
bounded_forward(b, theta)
```

## Arguments

- b:

  A length-2 numeric vector, `c(lower, upper)`.

- theta:

  A numeric vector strictly inside the bounds.

## Value

A numeric vector on the unconstrained scale, as long as `theta`. A value
on or outside a bound gives an infinite or `NaN` entry rather than an
error, the check belonging to
[`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md).

## Details

Strictly inside. A value **on** a bound maps to an infinite \\\eta\\, so
a run started there begins at infinity and fails far from its cause.
That is why
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
rejects such a starting value by name, before anything is evaluated.

## See also

[`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)
for the inverse and its derivatives,
[`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
for the validation.

## Examples

``` r
bounded_forward(c(0, Inf), c(0.5, 1, 8))
#> [1] -0.6931472  0.0000000  2.0794415
bounded_forward(c(0, 1), 0.5)
#> [1] 0

# A value on a bound has no finite image, so such a starting point is
# unusable and minimize() refuses it.
bounded_forward(c(0, 1), c(0, 1))
#> [1] -Inf  Inf

# Round trip, in this direction too.
eta <- c(-2, 0, 3)
all.equal(bounded_forward(c(0, Inf), bounded_transform(c(0, Inf), eta)$h),
          eta)
#> [1] TRUE
```
