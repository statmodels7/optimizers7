# Eigenvalue Modification of the Hessian

Builds the modification that decomposes the Hessian, \\H = Q \Lambda
Q^\top\\, and raises every eigenvalue below a floor to the floor:
\$\$\tilde H = Q \max(\Lambda, \lambda\_{\min}) Q^\top, \qquad
\lambda\_{\min} = \max(\texttt{floor},\\
\texttt{floor}\\\max_i\lvert\lambda_i\rvert).\$\$ The eigenvectors and
the eigenvalues above the floor are unchanged. The step component along
an eigenvector \\q_i\\ is \\-q_i^\top g / \max(\lambda_i,
\lambda\_{\min})\\: the Newton component where \\\lambda_i \ge
\lambda\_{\min}\\, and, where the eigenvalue was raised, a component of
length \\\lvert q_i^\top g\rvert/\lambda\_{\min}\\ that is large when
the floor is small. The cost is a symmetric eigendecomposition. It is
the default of
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md).
If the eigendecomposition fails, the ridge modification is used.

## Usage

``` r
modify_eigen(floor = 1e-08)
```

## Arguments

- floor:

  The lower bound of the eigenvalues, a single positive number. Defaults
  to `1e-8`. The bound is relative to the largest eigenvalue in absolute
  value when that value exceeds one.

## Value

An S7 object of class `ModifyEigen`, inheriting from
[`hessian_modification()`](https://statmodels7.github.io/optimizers7/reference/hessian_modification.md),
with the properties `floor` and `label`.

## See also

[`modify_ridge()`](https://statmodels7.github.io/optimizers7/reference/modify_ridge.md),
[modifications](https://statmodels7.github.io/optimizers7/reference/modifications.md)
for the comparison,
[`hessian_modify()`](https://statmodels7.github.io/optimizers7/reference/hessian_modify.md),
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)

## Examples

``` r
m <- modify_eigen(floor = 1e-6)
m
#> <hessian_modification> eigenvalue floor (floor = 1e-06)

# The method for hessian_modify() applied to an indefinite matrix.
H <- matrix(c(2, 1, 1, -3), 2, 2)
eigen(H, only.values = TRUE)$values
#> [1]  2.192582 -3.192582
eigen(hessian_modify(m, H), only.values = TRUE)$values
#> [1] 2.192582e+00 3.192582e-06

# Passed to the optimizer.
f <- function(p) p[1]^4 - 2 * p[1]^2 + p[2]^2
g <- function(p) c(4 * p[1]^3 - 4 * p[1], 2 * p[2])
h <- function(p) matrix(c(12 * p[1]^2 - 4, 0, 0, 2), 2, 2)
minimize(newton(modify = modify_eigen(floor = 1e-6)), f, c(0.05, 0.6),
         gr = g, hess = h)
#> <optimizer_result> Newton
#>   value      : -1
#>   par        : 1 0
#>   iterations : 4   evaluations: f 5, g 5
#>   elapsed    : 0 us
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
```
