# Nonlinear Conjugate Gradients

Gradient descent with the zigzag taken out: each direction is bent back
towards the previous one. It stores two vectors and no matrix, which
makes it the classical answer for a problem too large to hold a Hessian.

## Usage

``` r
cg(
  criterion = crit_any(crit_grad(1e-08), crit_rel_obj(1e-12)),
  beta = c("pr", "fr", "hs", "dy"),
  restart_every = 0,
  step = 1,
  line_search = wolfe(c2 = 0.1),
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

- beta:

  Which formula for the bend: `"pr"` (Polak–Ribière, the default),
  `"fr"` (Fletcher–Reeves), `"hs"` (Hestenes–Stiefel) or `"dy"`
  (Dai–Yuan). See Details.

- restart_every:

  Restart at steepest descent every this many iterations. Defaults to
  `0`, meaning never; a positive value is usually the dimension of the
  problem.

- step:

  Initial step length offered to the line search. Defaults to `1`.

- line_search:

  Defaults to `wolfe(c2 = 0.1)`, and see Details — neither the search
  nor that constant is a free choice.

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

An S7 object of class `Cg`, inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Details

The direction is \$\$d_k = -g_k + \beta_k d\_{k-1},\$\$ and the whole
method is the choice of \\\beta\\. On a quadratic with an exact line
search the directions come out conjugate with respect to the Hessian, so
the method terminates in \\p\\ steps exactly — *without ever forming
that Hessian*, which is the point. The storage is two vectors against
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)'s
\\p \times p\\ matrix.

### Which beta

All four agree on a quadratic with an exact line search and differ
everywhere else. `"fr"` has the cleanest convergence theory and the
well-known practical fault of stalling for many iterations after a poor
step. `"pr"` recovers from a poor step immediately, because a small \\y
= g_k - g\_{k-1}\\ sends \\\beta\\ towards zero and the method back to
steepest descent; it can fail to converge in theory, and the standard
repair is to clamp \\\beta\\ at zero, which is what this does. The clamp
is not a fudge — \\\beta = 0\\ *is* a restart, and it is reported as one
in the trace. `"hs"` and `"dy"` are the other two standard choices.

### The line search is not optional, and neither is its tolerance

The theory behind every one of these formulas assumes a step satisfying
the **strong** Wolfe conditions, and uses it to prove that the direction
produced is a descent direction at all. Backtracking gives no such
guarantee, so
[`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
is the default and departing from it is departing from the theory.

The constant matters as much as the search.
[`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
defaults to \\c_2 = 0.9\\, which is right for
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
where the curvature approximation repairs a loose step. Conjugate
gradients has nothing to repair with: the accumulated conjugacy is only
as good as the line search that produced it, and a loose one degrades it
directly. Measured on Rosenbrock, \\c_2 = 0.9\\ needs 120 iterations
against 35 at \\c_2 = 0.1\\, which is why the default here is the
tighter one.

As a safeguard against the cases the theory misses, a direction that
comes out non-descent is replaced by \\-g\\ and the substitution is
reported.

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

[`gd`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md),
[`bb`](https://statmodels7.github.io/optimizers7/reference/bb.md)

## Examples

``` r
cg()
#> <optimizer> conjugate gradients (pr)
#>   stop when : gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
#>   budgets   : maxit 1000, evaluations 20000
#>   settings  : beta = pr, restart_every = 0, step = 1, line_search = strong Wolfe (c1 = 1e-04, c2 = 0.1)
cg(beta = "fr", restart_every = 10)
#> <optimizer> conjugate gradients (fr)
#>   stop when : gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
#>   budgets   : maxit 1000, evaluations 20000
#>   settings  : beta = fr, restart_every = 10, step = 1, line_search = strong Wolfe (c1 = 1e-04, c2 = 0.1)

f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))
minimize(cg(), f, c(-1.2, 1), gr = gr)@par
#> [1] 1 1
```
