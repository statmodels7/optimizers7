# The Barzilai-Borwein Method

Gradient descent with one number changed, and the number is what makes
it work: the step length is a curvature estimate taken from the previous
step rather than found by searching.

## Usage

``` r
bb(
  criterion = crit_any(crit_grad(1e-08), crit_rel_obj(1e-12)),
  variant = c("alternate", "bb1", "bb2"),
  alpha0 = 0.01,
  alpha_min = 1e-10,
  alpha_max = 1e+10,
  step = 1,
  line_search = armijo(),
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

- step:

  Initial multiplier offered to the line search. Defaults to `1`, so the
  Barzilai-Borwein step is tried unaltered first.

- line_search:

  Defaults to
  [`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md);
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

It reaches every smooth minimum in
[`test_problems`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
to machine precision while storing two vectors and doing no linear
algebra at all, and it solves a quadratic in two iterations, because the
first secant pair already contains the answer. For a method with one
scalar of memory that is a great deal.

What it does not do is get there quickly on a curved valley. Measured on
Rosenbrock it takes about 930 iterations where
[`cg`](https://statmodels7.github.io/optimizers7/reference/cg.md) takes
35 and
[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/Lbfgs.md)
takes 35 — a thousandth of the memory and thirty times the iterations.
See below for why, since the reason is this implementation rather than
the method.

### Which variant

`"bb1"` and `"bb2"` are the two quotients, and `"alternate"` switches
between them at each iteration, which is the usual recommendation and
the default: they estimate the same curvature from different ends, and
taking them in turn is more robust than trusting either.

### The line search, and what it is for here

The step length is the method, so the line search must not be allowed to
reshape it. It is offered the Barzilai-Borwein step first and unaltered,
and only backtracks when that step fails the sufficient-decrease test.

That backtracking is where the iteration count above comes from, and it
is worth being plain about. The method is **famously non-monotone**: its
efficiency comes from steps that make the objective worse now in order
to align with the curvature, and an Armijo condition forbids exactly
those. The classical remedy is a *nonmonotone* line search, which
requires sufficient decrease against the worst of the last several
values rather than the last one, and with it Barzilai-Borwein is
competitive with conjugate gradients. This package has no such search,
so what is here is the safe version and substantially the slower one on
a curved valley. It is not a defect in the method and it is a real
limitation of this implementation.

A secant pair reporting non-positive curvature is skipped rather than
used, and a step length outside `alpha_min` and `alpha_max` is clamped.
Both are reported in the trace: an unbounded step is how a first-order
method leaves the domain entirely.

## References

Barzilai, J. and Borwein, J. M. (1988). Two-point step size gradient
methods. *IMA Journal of Numerical Analysis* **8**, 141–148.

## See also

[`gd`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`cg`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/Lbfgs.md)

## Examples

``` r
bb()
#> <optimizer> barzilai-borwein (alternate)
#>   stop when : gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
#>   budgets   : maxit 1000, evaluations 20000
#>   settings  : variant = alternate, alpha0 = 0.01, alpha_min = 1e-10, alpha_max = 1e+10, step = 1, line_search = Armijo backtracking (c1 = 1e-04)

f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))
minimize(bb(), f, c(-1.2, 1), gr = gr)@par
#> [1] 1 1

# two iterations on a quadratic, because one secant pair is the whole answer
minimize(bb(), function(p) sum((p - c(1, 2))^2), c(0, 0),
         gr = function(p) 2 * (p - c(1, 2)))@iterations
#> [1] 2
```
