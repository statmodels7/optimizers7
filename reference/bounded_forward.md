# The Forward Bound Transform

Maps a parameter from inside its bounds to the unconstrained scale: the
inverse of
[`bounded_transform`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)'s
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

A numeric vector on the unconstrained scale.

## Details

Strictly inside. A value on a bound maps to an infinite \\\eta\\, so a
run started there begins at infinity and fails far from its cause; that
is why
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
refuses such a starting value by name.

## See also

[`bounded_transform`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)

## Examples

``` r
bounded_forward(c(0, Inf), c(0.5, 1, 8))
#> [1] -0.6931472  0.0000000  2.0794415
bounded_forward(c(0, 1), 0.5)
#> [1] 0
```
