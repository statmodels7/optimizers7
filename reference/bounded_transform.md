# The Bound Transform, and Its First Two Derivatives

Evaluates the reparametrization defined by a pair of bounds. The package
removes a box in this way, and the function is exported so that a
user-defined
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
method can remove one in the same way.

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

Bounds are removed here and not enforced. The map \\h\\ is a shifted
exponential for a one-sided bound, a scaled logistic function for a
two-sided one, and the identity for neither (its inverse,
[`bounded_forward()`](https://statmodels7.github.io/optimizers7/reference/bounded_forward.md),
is the shifted log or the scaled logit), so optimizing in \\\eta\\ makes
every proposed point admissible by construction.

The image of \\\eta = 0\\ is the middle of the box in the sense of the
map: `1` for \\(0, \infty)\\, `0.5` for \\(0, 1)\\, `4` for \\(2, 6)\\,
and `0` where there is no bound. For this reason
[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
gives a sensible start for every kind of parameter.

## Using it in a method written by the caller

Map the starting value with
[`bounded_forward()`](https://statmodels7.github.io/optimizers7/reference/bounded_forward.md),
run unconstrained, and wrap the objective so that it maps back before
evaluating. The Jacobian is diagonal, so the chain rule is one product
per coordinate, \\\partial f/\partial \eta_i = (\partial f/\partial
\theta_i)\\ h_i'\\, and `d2` is needed only when a Hessian is
transformed, where it appears on the diagonal alone. Report `par` on the
caller's scale.

## Relation to linkfunctions7

\\h\\ and its two derivatives are the inverse link of `bounded_link()`
of linkfunctions7 and its derivatives, written in C++ because the
transform is applied at every evaluation of the objective, where a
callback into R would cost more than the loop saves. The test suite
compares them with `linkinv()`, `dlinkinv()` and `d2linkinv()`, so that
the copy stays equal to the original.

## See also

[`bounded_forward()`](https://statmodels7.github.io/optimizers7/reference/bounded_forward.md)
for the inverse,
[`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
for the shape in which `b` is supplied,
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
for the run that uses both.

## Examples

``` r
# For a variance, the whole line maps onto the positive half-line, and every
# value of eta, however large, maps inside it.
bounded_transform(c(0, Inf), c(-2, 0, 2))$h
#> [1] 0.1353353 1.0000000 7.3890561
bounded_transform(c(0, Inf), c(-500, 500))$h > 0
#> [1] TRUE TRUE

# For a probability, the result also includes the derivative that carries a
# gradient across.
str(bounded_transform(c(0, 1), c(-1, 0, 1)))
#> List of 3
#>  $ h : num [1:3] 0.269 0.5 0.731
#>  $ d1: num [1:3] 0.197 0.25 0.197
#>  $ d2: num [1:3] 0.0909 0 -0.0909

# Zero maps to the middle of the box, whatever the box. For this reason
# start_zeros() gives a sensible start for every kind of parameter.
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
