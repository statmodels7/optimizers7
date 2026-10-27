# Ridge Modification of the Hessian

Builds the modification that adds a multiple of the identity, \\\tilde H
= H + \tau I\\. The shift starts at \\\tau = \max(\texttt{floor},\\
\texttt{floor} - \min_i H\_{ii})\\, which makes the diagonal positive,
and is doubled until the Cholesky factorization of \\\tilde H\\ succeeds
(at most 60 times; the identity matrix is used if none succeeds). This
is the idea of Levenberg: the direction moves from the Newton direction
at \\\tau = 0\\ to a scaled steepest-descent direction for large
\\\tau\\. It costs less than
[`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md)
and shifts every eigenvalue, not only the small ones.

## Usage

``` r
modify_ridge(floor = 1e-08)
```

## Arguments

- floor:

  The smallest shift, a single positive number. Defaults to `1e-8`.

## Value

An S7 object of class `ModifyRidge`, inheriting from
[`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md),
with the properties `floor` and `label`.

## See also

[`modify_eigen()`](https://statmodels7.github.io/optimizers7/reference/modify_eigen.md),
[modifications](https://statmodels7.github.io/optimizers7/reference/modifications.md)
for the comparison,
[`hessian_modify()`](https://statmodels7.github.io/optimizers7/reference/hessian_modify.md),
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)

## Examples

``` r
m <- modify_ridge(floor = 1e-6)
m
#> <hessian_modification> ridge (floor = 1e-06)

# The method for hessian_modify() applied to an indefinite matrix.
H <- matrix(c(2, 1, 1, -3), 2, 2)
eigen(hessian_modify(m, H), only.values = TRUE)$values
#> [1] 8.192584 2.807420

# Passed to the optimizer.
f <- function(p) p[1]^4 - 2 * p[1]^2 + p[2]^2
g <- function(p) c(4 * p[1]^3 - 4 * p[1], 2 * p[2])
h <- function(p) matrix(c(12 * p[1]^2 - 4, 0, 0, 2), 2, 2)
minimize(newton(modify = modify_ridge(floor = 1e-6)), f, c(0.05, 0.6),
         gr = g, he = h)
#> <optimizer_result> Newton
#>   value      : -1
#>   par        : 1 0
#>   iterations : 4   evaluations: f 5, g 5
#>   elapsed    : 0 us
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
```
