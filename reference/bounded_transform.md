# The Bound Transform, and Its First Two Derivatives

Evaluates the reparametrisation a set of bounds implies. This is how the
package removes a box, and it is exported so that a
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
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

A list with `h` (the parameter), `d1` and `d2`.

## Details

Bounds are not enforced here, they are removed: a shifted log for a
one-sided bound, a scaled logit for two, and the identity for neither.
Optimise in \\\eta\\ and every proposed point is admissible by
construction.

To honour `bounds` in a user-defined method: map the starting value with
[`bounded_forward`](https://statmodels7.github.io/optimizers7/reference/bounded_forward.md),
run unconstrained, and wrap the objective so that it maps back before
evaluating. The Jacobian is diagonal, so the chain rule is short —
\\\partial f/\partial \eta_i = (\partial f/\partial \theta_i)\\ h_i'\\ —
and `d2` is needed only when a Hessian is transformed, where it appears
on the diagonal alone. Report `par` on the user's scale.

These are linkfunctions7's `bounded_link()`, written out in C++ because
the transform is applied on every objective evaluation and a callback
into R there would undo the reason for compiling the loop. The test
suite pins them to `linkinv()`, `dlinkinv()` and `d2linkinv()` on every
run, so the copy cannot drift from the original.

## See also

[`bounded_forward`](https://statmodels7.github.io/optimizers7/reference/bounded_forward.md),
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

## Examples

``` r
# a variance: the whole line maps onto the positive half
bounded_transform(c(0, Inf), c(-2, 0, 2))$h
#> [1] 0.1353353 1.0000000 7.3890561

# a probability, and the derivative that carries a gradient across
str(bounded_transform(c(0, 1), c(-1, 0, 1)))
#> List of 3
#>  $ h : num [1:3] 0.269 0.5 0.731
#>  $ d1: num [1:3] 0.197 0.25 0.197
#>  $ d2: num [1:3] 0.0909 0 -0.0909

# the round trip
bounded_transform(c(0, 1), bounded_forward(c(0, 1), c(0.1, 0.5, 0.9)))$h
#> [1] 0.1 0.5 0.9
```
