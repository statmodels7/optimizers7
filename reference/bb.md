# The Barzilai-Borwein Method

Gradient descent with the Barzilai-Borwein step length: the direction is
\\-\alpha_k g_k\\, where \\\alpha_k\\ is a scalar estimate of the
inverse curvature computed from the previous secant pair.

## Usage

``` r
bb(
  criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
  variant = c("alternate", "bb1", "bb2"),
  alpha0 = 0.01,
  alpha_min = 1e-10,
  alpha_max = 1e+10,
  curv_tol = 1e-10,
  step = 1,
  line_search = nonmonotone(),
  maxit = 1000,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 20,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule; see
  [`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md).

- variant:

  `"alternate"` (default), `"bb1"` or `"bb2"`; see Details.

- alpha0:

  The step length used on the first iteration, before there is a secant
  pair to estimate one from. Defaults to `1e-2`.

- alpha_min, alpha_max:

  Bounds on the step length. Defaults `1e-10` and `1e10`.

- curv_tol:

  The relative curvature threshold: a secant pair is rejected when
  \\s^\top y \le c \lVert s \rVert \lVert y \rVert\\ for \\c\\ equal to
  `curv_tol`. Defaults to `1e-10`, which is the same relative test
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  applies to the same quantity.

- step:

  Initial multiplier offered to the line search. Defaults to `1`, so the
  Barzilai-Borwein step is tried unaltered first.

- line_search:

  The acceptance test for a trial step. Defaults to
  [`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md);
  see Details.

- maxit:

  Maximum iterations. Defaults to 1000.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`: no evaluation
  budget, so the run stops on the criterion or on `maxit`. Set a finite
  value to cap the cost of a run.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 20.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

## Value

An S7 object of class
[Bb](https://statmodels7.github.io/optimizers7/reference/Bb-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Details

Take the direction \\-\alpha g\\ with \$\$\alpha\_{BB1} = \frac{s^\top
s}{s^\top y}, \qquad \alpha\_{BB2} = \frac{s^\top y}{y^\top y},\$\$ the
two Rayleigh quotients of the secant pair \\s = x_k - x\_{k-1}\\, \\y =
g_k - g\_{k-1}\\. Both estimate the inverse curvature along the
direction just traveled, so this is a quasi-Newton method that has
discarded everything except one scalar. On a quadratic, where the
curvature is constant, that scalar is exactly right.

On a quadratic, where the curvature is constant, the estimate is exact
and the method converges in **two** iterations. On a general smooth
objective it needs more iterations than
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
while storing a single scalar instead of a matrix: on Rosenbrock 58
iterations against 35, but 67 objective evaluations against 49.

## Variants

`"bb1"` and `"bb2"` are the two quotients above, and `"alternate"`, the
default, switches between them at each iteration. They estimate the same
curvature from opposite ends, and alternating is the more robust choice
across problems, though not always the fastest on any one: measured on
Rosenbrock, `bb1` takes 56 iterations and 97 evaluations, `bb2` 51 and
56, `alternate` 58 and 67.

## Line search

The Barzilai-Borwein step is offered to the line search first and
unaltered, and backtracking occurs only when it fails the acceptance
test. The method makes progress through steps that may increase the
objective temporarily, so the default acceptance test is
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md),
which asks for improvement over the maximum of the last `memory` values
instead of over the current one. A plain Armijo condition rejects
exactly the steps the method relies on: measured on Rosenbrock, 72
iterations and 154 evaluations against 58 and 67.

`nonmonotone(memory = 0)` is
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
value for value, and the two give the identical run here, 72 iterations
and 154 evaluations. That identity is what makes the comparison a
comparison of the memory alone.

## Rejected secant pairs

A pair is used only if it reports positive curvature by a relative
margin, \\s^\top y \> c \lVert s \rVert \lVert y \rVert\\ with \\c\\ the
`curv_tol` argument, the same test
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
applies. When a pair is rejected the step length is reset to \\1/\lVert
g \rVert\_\infty\\, a trial displacement of order one in the parameters.
The reset reads the current gradient and not the step length being
replaced or a fixed constant, so it can neither freeze the iteration at
a too-short step nor produce one the backtracking cannot rescale. Steps
outside `[alpha_min, alpha_max]` are clamped, and both the reset and the
clamp appear in the trace, as `bb curvature reset` and `step shortened`.

## References

Barzilai, J. and Borwein, J. M. (1988). Two-point step size gradient
methods. *IMA Journal of Numerical Analysis* **8**, 141–148.

## See also

[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) for
the same direction with a line-searched step,
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
for the acceptance test this method needs,
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
for the next amount of curvature to carry.

## Examples

``` r
bb()
#> <optimizer> barzilai-borwein (alternate)
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 1000, evaluations Inf
#>   settings  : variant = alternate, alpha0 = 0.01, alpha_min = 1e-10, alpha_max = 1e+10, curv_tol = 1e-10, step = 1, line_search = nonmonotone backtracking (memory = 10)

f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))
minimize(bb(), f, c(-1.2, 1), gr = gr)@par
#> [1] 0.9999277 0.9998551

# Two iterations on a quadratic: one secant pair determines the curvature,
# and on a quadratic that curvature is exactly right.
minimize(bb(), function(p) sum((p - c(1, 2))^2), c(0, 0),
         gr = function(p) 2 * (p - c(1, 2)))@iterations
#> [1] 2

# The nonmonotone rule is what the method needs. With memory = 0 it is
# armijo() value for value, and both give the identical, slower run.
evals <- function(ls) minimize(bb(line_search = ls), f, c(-1.2, 1),
                               gr = gr)@counts[["f"]]
c(nonmonotone = evals(nonmonotone()),
  memory_zero = evals(nonmonotone(memory = 0)),
  armijo      = evals(armijo()))
#> nonmonotone memory_zero      armijo 
#>          67         154         154 

# Which variant is fastest is a property of the problem.
vapply(c("bb1", "bb2", "alternate"),
       function(v) minimize(bb(variant = v), f, c(-1.2, 1), gr = gr)@iterations,
       integer(1))
#>       bb1       bb2 alternate 
#>        56        51        58 
```
