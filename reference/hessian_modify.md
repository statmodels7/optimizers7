# Modify a Hessian

Generic that returns a positive definite matrix \\\tilde H\\ in place of
a Hessian \\H\\ for which the Cholesky factorization failed.
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
solve \\\tilde H d = -g\\ for the direction. There is one method for
each
[`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md)
subclass; a new modification is a subclass and a method, and
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
needs no change.

## Usage

``` r
hessian_modify(modify, H)
```

## Arguments

- modify:

  A
  [`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md)
  object.

- H:

  A symmetric numeric matrix.

## Value

A numeric matrix with the dimension of `H`, intended to be positive
definite.

## Requirements on a method

The method receives the modification object and the Hessian \\H\\, a
symmetric numeric matrix with finite entries, and returns a numeric
matrix of the same dimension. The returned matrix should be positive
definite: the method of the optimizer solves with it, and a direction
that is not a descent direction is replaced by a scaled gradient step,
which is recorded in the trace as `hessian solve failed`. A returned
matrix of another dimension or with non-finite entries is treated in the
same way. The method is called only at an iterate at which \\H\\ is not
positive definite, and is called from compiled code through R, once per
such iterate.

## Writing a modification

A modification is defined by (1) a subclass of
[`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md)
with the parameters of the rule as properties, (2) a constructor that
sets the `label` shown when the optimizer is printed, and (3) a method
for this generic, registered with
[`S7::method()`](https://rconsortium.github.io/S7/reference/method.html).
The shipped classes
[`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md)
and
[`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md)
are handled in compiled code; any other class, including a subclass of
one of them, is called through this generic.

## References

Dauphin, Y. N., Pascanu, R., Gulcehre, C., Cho, K., Ganguli, S. and
Bengio, Y. (2014). Identifying and attacking the saddle point problem in
high-dimensional non-convex optimization. *Advances in Neural
Information Processing Systems* **27**.

## See also

[`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md),
[`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md),
[`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md),
[modifications](https://statmodels7.github.io/optimizers7/reference/modifications.md),
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)

## Examples

``` r
# The shipped methods, on an indefinite matrix with eigenvalues 2 and -1.
H <- matrix(c(0.5, 1.5, 1.5, 0.5), 2, 2)
eigen(H, only.values = TRUE)$values
#> [1]  2 -1
hessian_modify(modify_eigen(floor = 0.1), H)
#>      [,1] [,2]
#> [1,]  1.1  0.9
#> [2,]  0.9  1.1
hessian_modify(modify_ridge(), H)
#>          [,1]     [,2]
#> [1,] 1.842177 1.500000
#> [2,] 1.500000 1.842177

# A user-defined modification: replace each eigenvalue by its absolute
# value, with a lower bound (the saddle-free Newton method of Dauphin et
# al., 2014). Unlike the eigenvalue floor, it uses the magnitude of a
# negative curvature as the curvature of the quadratic model.
ModifyAbs <- S7::new_class("ModifyAbs", parent = hessian_modification,
  properties = list(floor = S7::class_numeric))

modify_abs <- function(floor = 1e-8) {
  ModifyAbs(label = paste0("absolute eigenvalues (floor = ", floor, ")"),
            floor = floor)
}

S7::method(hessian_modify, ModifyAbs) <- function(modify, H) {
  e <- eigen(H, symmetric = TRUE)
  e$vectors %*% (pmax(abs(e$values), modify@floor) * t(e$vectors))
}

hessian_modify(modify_abs(), H)
#>      [,1] [,2]
#> [1,]  1.5  0.5
#> [2,]  0.5  1.5

# The new object is passed to newton() like a shipped one, here on the
# Himmelblau function, started where its Hessian is not positive definite.
k <- test_problems("himmelblau")[[1]]
eigen(k$he(c(0.5, 0.5)), only.values = TRUE)$values
#> [1] -20.05573 -37.94427
r <- minimize(newton(modify = modify_abs(), keep_trace = TRUE), k$fn,
              c(0.5, 0.5), gr = k$gr, he = k$he)
c(value = r@value, iterations = r@iterations, converged = r@converged)
#>        value   iterations    converged 
#> 3.346705e-24 7.000000e+00 1.000000e+00 
unique(r@trace$safeguard)
#> [1] "hessian modified (capped)" "none"                     
```
