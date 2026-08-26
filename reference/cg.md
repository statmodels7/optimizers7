# Nonlinear Conjugate Gradients

Nonlinear conjugate gradients: the direction is \\d_k = -g_k + \beta_k
d\_{k-1}\\, so consecutive directions are conjugate on a quadratic
rather than orthogonal. Storage is two vectors, with no matrix, which
suits problems too large for a Hessian.

## Usage

``` r
cg(
  criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
  beta = c("pr", "fr", "hs", "dy"),
  restart_every = 0,
  step = 1,
  line_search = wolfe(c2 = 0.1),
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

- beta:

  Which formula for the bend: `"pr"` (Polak-Ribière, the default),
  `"fr"` (Fletcher-Reeves), `"hs"` (Hestenes-Stiefel) or `"dy"`
  (Dai-Yuan). See Details.

- restart_every:

  Restart at steepest descent every this many iterations. Defaults to
  `0`, meaning never; a positive value is usually the dimension of the
  problem.

- step:

  Initial step length offered to the line search. Defaults to `1`.

- line_search:

  Defaults to `wolfe(c2 = 0.1)`; the convergence theory of the method
  assumes a strong Wolfe step, see Details.

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
[Cg](https://statmodels7.github.io/optimizers7/reference/Cg-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Details

The direction is \$\$d_k = -g_k + \beta_k d\_{k-1},\$\$ and the methods
differ only in the choice of \\\beta\\. The storage is two vectors
against
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)'s
\\p \times p\\ matrix.

## What conjugacy buys, and what it needs

On a quadratic **with an exact line search** the directions come out
conjugate with respect to the Hessian, and the method then terminates in
\\p\\ steps without ever forming that Hessian. The line search here is
inexact, so that guarantee does not transfer: measured on dense
quadratics at `crit_grad(1e-12)`, the run takes 16, 34 and 70 iterations
at \\p = 3, 8, 20\\. The saving is still real, since each iteration
costs two vectors.

## Choice of beta

All four agree on a quadratic under an exact line search and differ
everywhere else. `"fr"` has the cleanest convergence theory and the
well-known practical fault of stalling for many iterations after a poor
step. `"pr"` recovers from a poor step immediately, a small \\y = g_k -
g\_{k-1}\\ sending \\\beta\\ towards zero and the method back to
steepest descent; its known theoretical non-convergence is repaired by
clamping \\\beta\\ at zero, which restarts the method and appears in the
trace as `cg restart`. `"hs"` and `"dy"` are the other two standard
choices.

Measured on Rosenbrock from the customary start, iterations and
objective evaluations: `pr` 25 and 249 with twelve clamps, `fr` 59 and
743, `hs` 17 and 166, `dy` 50 and 667. `pr` is the default as the safest
of the four on a general objective, and `hs` was the fastest here; on
another problem the order will differ.

## Line search

The theory behind every one of these formulas assumes a step satisfying
the **strong** Wolfe conditions, and uses it to prove that the direction
produced is a descent direction at all. Backtracking gives no such
guarantee, so
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
is the default and departing from it departs from the theory.

The constant matters too.
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
defaults to \\c_2 = 0.9\\, which is right for
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
where the curvature approximation repairs a loose step. Conjugate
gradients has nothing to repair with: the accumulated conjugacy is only
as good as the line search that produced it. The default here is
therefore `wolfe(c2 = 0.1)`, and the cost of loosening it is modest on
Rosenbrock, 28 iterations against 25, with the value reached `3.3e-08`
against `7.4e-10`.

A direction that comes out non-descent is replaced by \\-g\\ and the
substitution is reported in the trace, as a safeguard against the cases
the theory misses.

## References

Hestenes, M. R. and Stiefel, E. (1952). Methods of conjugate gradients
for solving linear systems. *Journal of Research of the NBS* **49**,
409–436.

Polak, E. and Ribière, G. (1969). Note sur la convergence de méthodes de
directions conjuguées. *Revue Française d'Informatique et de Recherche
Opérationnelle* **3**, 35–43.

Dai, Y. H. and Yuan, Y. (1999). A nonlinear conjugate gradient method
with a strong global convergence property. *SIAM Journal on
Optimization* **10**, 177–182.

## See also

[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) for
the direction this one bends,
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
for a method with the same storage order and more curvature,
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) for
the scalar estimate.

## Examples

``` r
cg()
#> <optimizer> conjugate gradients (pr)
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 1000, evaluations Inf
#>   settings  : beta = pr, restart_every = 0, step = 1, line_search = strong Wolfe (c1 = 1e-04, c2 = 0.1)
cg(beta = "fr", restart_every = 10)
#> <optimizer> conjugate gradients (fr)
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 1000, evaluations Inf
#>   settings  : beta = fr, restart_every = 10, step = 1, line_search = strong Wolfe (c1 = 1e-04, c2 = 0.1)

f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))
minimize(cg(), f, c(-1.2, 1), gr = gr)@par
#> [1] 0.9999728 0.9999455

# The four formulas on the same problem. They agree on a quadratic under an
# exact line search and differ here, and which is best is a property of the
# problem.
vapply(c("pr", "fr", "hs", "dy"),
       function(b) minimize(cg(beta = b), f, c(-1.2, 1), gr = gr)@iterations,
       integer(1))
#> pr fr hs dy 
#> 25 59 17 50 

# Polak-Ribiere clamps beta at zero when a step goes badly, which restarts
# the method; the trace counts those.
r <- minimize(cg(keep_trace = TRUE), f, c(-1.2, 1), gr = gr)
table(r@trace$safeguard)
#> 
#>    cg restart step adjusted 
#>            12            13 
```
