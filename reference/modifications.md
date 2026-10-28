# Modifications of the Hessian

Overview of the objects that repair a Hessian that is not positive
definite, passed to
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
as `modify`.

## Value

Not applicable: this page documents the family. Each constructor returns
an object inheriting from
[`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md).

## Details

The Newton direction \\d = -H^{-1}g\\ is a descent direction only when
\\H\\ is positive definite. When the Cholesky factorization of \\H\\
fails, the method replaces \\H\\ by a positive definite matrix \\\tilde
H\\ and solves \\\tilde H d = -g\\. The object in `modify` selects how
\\\tilde H\\ is formed. The default is
[`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md).

- [`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md):

  Decomposes \\H = Q \Lambda Q^\top\\ and raises every eigenvalue below
  \\\max(\texttt{floor}, \texttt{floor}\\\max_i\lvert\lambda_i\rvert)\\
  to that value. It costs a symmetric eigendecomposition and leaves the
  eigenvectors and the eigenvalues above the floor unchanged. The step
  component along an eigenvector \\q_i\\ is \\-q_i^\top g /
  \max(\lambda_i, \lambda\_{\min})\\, the Newton component where
  \\\lambda_i \ge \lambda\_{\min}\\ and a component of length \\\lvert
  q_i^\top g\rvert/\lambda\_{\min}\\ where the eigenvalue was raised,
  which is large when the floor is small.

- [`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md):

  Adds \\\tau I\\ to \\H\\, starting from the smallest \\\tau \ge
  \texttt{floor}\\ that makes the diagonal positive and doubling
  \\\tau\\ until the Cholesky factorization succeeds. It costs a few
  factorizations and shifts every eigenvalue by the same amount.

- A user-defined subclass of
  [`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md):

  Any other rule, supplied as a method for
  [`hessian_modify()`](https://statmodels7.github.io/optimizers7/reference/hessian_modify.md).
  See that page for an executable example.

Both shipped modifications run in compiled code. A user-defined class is
called from the compiled loop through R at each iteration at which the
Hessian is not positive definite.

A repair is recorded in the trace as `hessian modified` (or
`hessian modified (capped)`) whichever modification was used; see
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

## See also

[`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md),
[`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md),
[`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md),
[`hessian_modify()`](https://statmodels7.github.io/optimizers7/reference/hessian_modify.md),
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)

## Examples

``` r
f <- function(p) p[1]^4 - 2 * p[1]^2 + p[2]^2
g <- function(p) c(4 * p[1]^3 - 4 * p[1], 2 * p[2])
h <- function(p) matrix(c(12 * p[1]^2 - 4, 0, 0, 2), 2, 2)

# The Hessian at the start has a negative eigenvalue.
eigen(h(c(0.05, 0.6)), only.values = TRUE)$values
#> [1]  2.00 -3.97

sapply(list(eigen = modify_eigen(), ridge = modify_ridge()), function(m) {
  r <- minimize(newton(modify = m), f, c(0.05, 0.6), gr = g, hess = h)
  c(value = r@value, iterations = r@iterations)
})
#>            eigen ridge
#> value         -1    -1
#> iterations     4     4
```
