# The Bound Transform, and Its First Two Derivatives

Evaluates the reparametrization a set of bounds implies. This is how the
package removes a box, and it is exported so that a
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
user-defined method can remove one the same way.

## Usage

``` r
bounded_transform(b, eta)
```

## Arguments

- b:

  A length-2 numeric vector, `c(lower, upper)`, using `-Inf` and `Inf`
  for a side that is unbounded.

- eta:

  A numeric vector on the unconstrained scale.

## Value

A named list of three numeric vectors, each as long as `eta`:

- `h`:

  the parameter \\\theta = h(\eta)\\, inside the box.

- `d1`:

  \\h'(\eta)\\, the chain-rule factor for a gradient.

- `d2`:

  \\h''(\eta)\\, needed only for a Hessian.

## Details

Bounds are removed here, not enforced. The map is a shifted log for a
one-sided bound, a scaled logit for a two-sided one, and the identity
for neither, so optimizing in \\\eta\\ makes every proposed point
admissible by construction.

What \\\eta = 0\\ becomes is the midpoint of the box in the map's own
sense: `1` for \\(0, \infty)\\, `0.5` for \\(0, 1)\\, `4` for \\(2,
6)\\, and `0` where there is no bound at all. That is why
[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
is a sensible default start whatever the parameter means.

## Using it in a method of your own

Map the starting value with
[`bounded_forward()`](https://statmodels7.github.io/optimizers7/reference/bounded_forward.md),
run unconstrained, and wrap the objective so that it maps back before
evaluating. The Jacobian is diagonal, so the chain rule is one product
per coordinate, \\\partial f/\partial \eta_i = (\partial f/\partial
\theta_i)\\ h_i'\\, and `d2` is needed only when a Hessian is
transformed, where it appears on the diagonal alone. Report `par` on the
caller's scale.

## Where it comes from

These are linkfunctions7's `bounded_link()` written out in C++, because
the transform is applied at every objective evaluation and a callback
into R there would undo the reason for compiling the loop. The test
suite pins them to `linkinv()`, `dlinkinv()` and `d2linkinv()` on every
run, so the copy cannot drift from the original.

## See also

[`bounded_forward()`](https://statmodels7.github.io/optimizers7/reference/bounded_forward.md)
for the inverse,
[`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
for the shape `b` comes in,
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
for the run that uses both.

## Examples

``` r
# A variance: the whole line maps onto the positive half, and no value of
# eta, however absurd, leaves it.
bounded_transform(c(0, Inf), c(-2, 0, 2))$h
#> [1] 0.1353353 1.0000000 7.3890561
bounded_transform(c(0, Inf), c(-500, 500))$h > 0
#> [1] TRUE TRUE

# A probability, with the derivative that carries a gradient across.
str(bounded_transform(c(0, 1), c(-1, 0, 1)))
#> List of 3
#>  $ h : num [1:3] 0.269 0.5 0.731
#>  $ d1: num [1:3] 0.197 0.25 0.197
#>  $ d2: num [1:3] 0.0909 0 -0.0909

# Zero is the middle of the box, whichever box it is. This is why
# start_zeros() means something sensible for every kind of parameter.
vapply(list(c(-Inf, Inf), c(0, Inf), c(0, 1), c(2, 6)),
       function(b) bounded_transform(b, 0)$h, 0)
#> [1] 0.0 1.0 0.5 4.0

# The two maps are inverse to each other.
theta <- c(0.1, 0.5, 0.9)
all.equal(bounded_transform(c(0, 1), bounded_forward(c(0, 1), theta))$h,
          theta)
#> [1] TRUE

# d1 is a derivative, and a central difference confirms it.
h <- 1e-6
(bounded_transform(c(0, 1), 0.3 + h)$h - bounded_transform(c(0, 1), 0.3 - h)$h) /
  (2 * h)
#> [1] 0.2444583
bounded_transform(c(0, 1), 0.3)$d1
#> [1] 0.2444583
```
