# The Barzilai-Borwein Method

Gradient descent whose step length is a scalar curvature estimate
computed from the previous step's secant pair rather than found by
searching.

## Usage

``` r
bb(
  criterion = crit_any(crit_grad(1e-08), crit_rel_obj(1e-12)),
  variant = c("alternate", "bb1", "bb2"),
  alpha0 = 0.01,
  alpha_min = 1e-10,
  alpha_max = 1e+10,
  curv_tol = 1e-10,
  step = 1,
  line_search = nonmonotone(),
  maxit = 1000,
  max_eval = 20000,
  verbose = FALSE,
  refresh = 20,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule; see
  [`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md).

- variant:

  `"alternate"` (default), `"bb1"` or `"bb2"`; see Details.

- alpha0:

  The step length used on the first iteration, before there is a secant
  pair to estimate one from. Defaults to `1e-2`.

- alpha_min, alpha_max:

  Bounds on the step length. Defaults `1e-10` and `1e10`.

- curv_tol:

  The relative curvature threshold: a secant pair is refused when
  \\s^\top y \le c \lVert s \rVert \lVert y \rVert\\ for \\c\\ equal to
  `curv_tol`. Defaults to `1e-10`, which is the same relative test
  [`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  applies to the same quantity.

- step:

  Initial multiplier offered to the line search. Defaults to `1`, so the
  Barzilai-Borwein step is tried unaltered first.

- line_search:

  The acceptance test for a trial step. Defaults to
  [`nonmonotone`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md);
  see Details.

- maxit:

  Maximum iterations. Defaults to 1000.

- max_eval:

  Maximum objective evaluations. Defaults to 20000.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 20.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

## Value

An S7 object of class `Bb`, inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Details

Take the direction \\-\alpha g\\ with \$\$\alpha\_{BB1} = \frac{s^\top
s}{s^\top y}, \qquad \alpha\_{BB2} = \frac{s^\top y}{y^\top y},\$\$ the
two Rayleigh quotients of the secant pair \\s = x_k - x\_{k-1}\\, \\y =
g_k - g\_{k-1}\\. Both estimate the inverse curvature along the
direction just travelled, so this is a quasi-Newton method that has
discarded everything except one scalar. On a quadratic, where the
curvature is constant, that scalar is exactly right.

On a quadratic, where the curvature is constant, the estimate is exact
and the method converges in two iterations; on a general smooth
objective it typically needs more iterations than
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
while storing a single scalar instead of a matrix.

### Variants

`"bb1"` and `"bb2"` are the two quotients above, and `"alternate"`, the
default, switches between them at each iteration: they estimate the same
curvature from opposite ends, and alternating them is more robust than
either alone.

### Line search

The Barzilai-Borwein step is offered to the line search first and
unaltered, and backtracking occurs only when it fails the acceptance
test. Because the method makes progress through steps that may increase
the objective temporarily, the default acceptance test is
[`nonmonotone`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md),
which requires improvement over the maximum of the last `memory` values
rather than over the current one; a plain Armijo condition rejects
exactly the steps the method relies on and slows it considerably.
[`nonmonotone`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
with `memory = 0` coincides with
[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md).

### Refused secant pairs

A pair is used only if it reports positive curvature by a relative
margin, \\s^\top y \> c \lVert s \rVert \lVert y \rVert\\ with \\c\\ the
`curv_tol` argument – the same test
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
applies. When a pair is refused, the step length is reset to \\1/\lVert
g \rVert\_\infty\\, giving a trial displacement of order one in the
parameters. This reset depends on the current gradient rather than on
the step length being replaced or on a fixed constant, so it can neither
freeze the iteration at a too-short step nor produce a step the
backtracking cannot rescale. Steps outside `[alpha_min, alpha_max]` are
clamped, and both the reset and the clamp are recorded in the trace.

## References

Barzilai, J. and Borwein, J. M. (1988). Two-point step size gradient
methods. *IMA Journal of Numerical Analysis* **8**, 141–148.

## See also

[`gd`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`cg`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)

## Examples

``` r
bb()
#> <optimizer> barzilai-borwein (alternate)
#>   stop when : gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
#>   budgets   : maxit 1000, evaluations 20000
#>   settings  : variant = alternate, alpha0 = 0.01, alpha_min = 1e-10, alpha_max = 1e+10, curv_tol = 1e-10, step = 1, line_search = nonmonotone backtracking (memory = 10)

f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))
minimize(bb(), f, c(-1.2, 1), gr = gr)@par
#> [1] 1 1

# two iterations on a quadratic: one secant pair determines the curvature
minimize(bb(), function(p) sum((p - c(1, 2))^2), c(0, 0),
         gr = function(p) 2 * (p - c(1, 2)))@iterations
#> [1] 2
```
