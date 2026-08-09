# Line Search Satisfying the Strong Wolfe Conditions

Sufficient decrease and a curvature condition together: \$\$f(x + s d)
\leq f(x) + c_1 s\\ g^\top d, \qquad \lvert g(x + s d)^\top d \rvert
\leq c_2 \lvert g^\top d \rvert .\$\$

## Usage

``` r
wolfe(c1 = 1e-04, c2 = 0.9, max_step = 30)
```

## Arguments

- c1:

  Sufficient-decrease constant. Defaults to `1e-4`.

- c2:

  Curvature constant, with \\c_1 \< c_2 \< 1\\. Defaults to `0.9`, the
  usual choice for a quasi-Newton method; `0.1` is conventional for
  nonlinear conjugate gradients, which want a more exact line search.

- max_step:

  Maximum trial steps in each of the bracketing and zoom phases.
  Defaults to 30.

## Value

A
[`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
object.

## Details

The curvature condition is what
[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
cannot provide, and it is not a refinement: a quasi-Newton method needs
it to work at all. BFGS builds its approximation from the secant pair
\\(s, y)\\ with \\y = g\_{new} - g\_{old}\\, and a step so short that
the gradient has barely moved gives a pair carrying no curvature
information — so the update must either be skipped or it corrupts the
matrix. Requiring the gradient along the direction to have shrunk by a
factor \\c_2\\ is exactly the guarantee that this does not happen.

The implementation is the bracketing-and-zoom scheme, with bisection
inside the zoom rather than polynomial interpolation: a few more
evaluations, and it cannot be defeated by an awkwardly shaped interval.

It costs gradient evaluations at trial points, which
[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
does not, so it is the more expensive choice per iteration and usually
the cheaper one per problem.

## References

Wolfe, P. (1969). Convergence conditions for ascent methods. *SIAM
Review* **11**, 226–235.

Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*, 2nd
edition. Springer, New York.

## See also

[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)

## Examples

``` r
wolfe()
#> <line_search> strong Wolfe (c1 = 1e-04, c2 = 0.9)
minimize(gd(line_search = wolfe(), maxit = 200),
         function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2, c(-1.2, 1))
#> <optimizer_result> gradient descent
#>   value      : 2.34837e-05
#>   par        : 1.0048 1.0097
#>   iterations : 200   evaluations: f 2796, g 0
#>   elapsed    : 16 ms
#>   converged  : NO (iteration budget reached)
#>   note       : gradient obtained by finite differences
```
