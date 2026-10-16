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

  Which formula for \\\beta\\: `"pr"` (Polak-Ribière, the default),
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
  assumes a Wolfe step, see Details.

- maxit:

  Maximum iterations. Defaults to 1000.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`, which sets no
  evaluation budget, so the run stops on the criterion or on `maxit`.
  Set a finite value to cap the cost of a run.

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

## Conjugate directions and the line search

On a quadratic with an exact line search the directions are conjugate
with respect to the Hessian, and the method terminates in \\p\\ steps
without forming that Hessian. The line search here is inexact, so that
guarantee does not hold, and the number of iterations on a quadratic
depends on its conditioning. Each iteration still stores only two
vectors.

## Choice of beta

All four formulas agree on a quadratic under an exact line search and
differ elsewhere. `"fr"` has the simplest convergence theory and can
stall for many iterations after a poor step. `"pr"` recovers from a poor
step at once, since a small \\y = g_k - g\_{k-1}\\ makes \\\beta\\ small
and the direction close to steepest descent; its known theoretical
non-convergence is removed by clamping \\\beta\\ at zero, which restarts
the method and appears in the trace as `cg restart`. `"hs"` and `"dy"`
are the other two standard choices. `"pr"` is the default as the safest
of the four on a general objective; which formula needs the fewest
iterations depends on the problem (see the examples).

## Line search

The convergence results for these formulas assume a step that satisfies
the Wolfe conditions (the strong ones for `"fr"` and `"pr"`) and use
them to show that the direction produced is a descent direction.
Backtracking gives no such guarantee, so
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
is the default.

[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
defaults to \\c_2 = 0.9\\, which suits
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
because its curvature update compensates for an inexact step. The
conjugacy accumulated by conjugate gradients depends directly on the
accuracy of the line search, so the default here is `wolfe(c2 = 0.1)`.

A direction that is not a descent direction is replaced by \\-g\\, and
the trace records the replacement as `cg not descent`. A periodic
restart set by `restart_every` is recorded as `cg restart`.

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
steepest descent,
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
for a method with storage of the same order and more curvature
information,
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

# The four formulas are compared on the same problem. They agree on a quadratic under
# an exact line search and differ here; which is best depends on the
# problem.
vapply(c("pr", "fr", "hs", "dy"),
       function(b) minimize(cg(beta = b), f, c(-1.2, 1), gr = gr)@iterations,
       integer(1))
#> pr fr hs dy 
#> 25 59 17 50 

# Polak-Ribiere clamps beta at zero after a poor step, which restarts the
# method; the trace counts the restarts.
r <- minimize(cg(keep_trace = TRUE), f, c(-1.2, 1), gr = gr)
table(r@trace$safeguard)
#> 
#>    cg restart step adjusted 
#>            12            13 
```
