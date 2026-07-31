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

  Defaults to
  [`nonmonotone`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md),
  and that is not an incidental choice; see Details.

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
first secant pair already contains the answer.

On Rosenbrock it takes 68 iterations and 77 objective evaluations,
against 65 for
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and
[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
and 353 for
[`cg`](https://statmodels7.github.io/optimizers7/reference/cg.md) — from
one stored number rather than a matrix, a list of secant pairs or a
direction. For a method of this size that is a remarkable place to be.

### Which variant

`"bb1"` and `"bb2"` are the two quotients, and `"alternate"` switches
between them at each iteration, which is the usual recommendation and
the default: they estimate the same curvature from different ends, and
taking them in turn is more robust than trusting either.

### The line search, and what it is for here

The step length is the method, so the line search must not be allowed to
reshape it. The Barzilai-Borwein step is offered first and unaltered,
and backtracking happens only when it fails the acceptance test.

Which test matters. The method is **non-monotone by design**: its
efficiency comes from steps that make the objective worse now in order
to be aligned with the curvature later, and an Armijo condition forbids
exactly those. The default is therefore
[`nonmonotone`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md),
which asks for improvement on the worst of the last several values
rather than on the present one. Measured on Rosenbrock it accepts eleven
uphill steps out of sixty-seven, where
[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
accepts none, and finishes in 68 iterations and 77 evaluations against
82 and 186. On the non-smooth problem in
[`test_problems`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
it ends eight orders of magnitude closer, though it takes 439 iterations
to get there against 29 that stop early.

[`nonmonotone`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
with `memory = 0` is exactly
[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
– the same trace, value for value – so the two can be compared without
changing anything else.

### What happens when a pair carries no curvature

A secant pair is usable only if it reports positive curvature, and the
test is stated relatively – \\s^\top y \le c \lVert s \rVert \lVert y
\rVert\\, with \\c\\ the `curv_tol` argument, the same test
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
applies to the same quantity – so that a value positive only by rounding
does not pass. When the pair is refused the step length is set to
\\1/\lVert g \rVert\_\infty\\, which asks for a trial displacement of
order one in the parameters whatever the gradient's magnitude.

That looks like an odd thing to reach for until one sees what the
obvious alternatives do. Keeping the previous \\\alpha\\ is an absorbing
state: a short step samples curvature over a short interval, on a valley
floor that curvature is negative, so the pair is refused, so \\\alpha\\
stays short and the next pair is refused too. Measured on Rosenbrock the
`"bb2"` variant fell into it at iteration six and never left: 873
refusals in 945 iterations, the gradient norm frozen while the objective
crept down by 0.002 a step. Restarting at `alpha0` escapes that trap and
then fails in the same shape wherever `alpha0` is itself too short to
escape with – on a boxed quadratic seen through its reparametrisation it
cost 1395 refusals in 1521 iterations, against 113 for keeping. Reaching
for `alpha_max`, as the spectral projected gradient literature does in a
setting that also projects and clamps differently, asks the line search
to backtrack a direction of length `1e10`, which thirty halvings do not
bring back into scale: on that same boxed problem the run stopped after
eight iterations a whole unit from the solution while reporting success.
Of the four, \\1/\lVert g \rVert\_\infty\\ is the only one that is never
the worst, over Rosenbrock, Beale, Powell's quartic and the boxed
quadratic; it cannot freeze, since it does not depend on the \\\alpha\\
it replaces, and it cannot explode, since it scales with the gradient.

A step length outside `alpha_min` and `alpha_max` is clamped. Both the
reset and the clamp are named in the trace: an unbounded step is how a
first-order method leaves the domain entirely.

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

# two iterations on a quadratic, because one secant pair is the whole answer
minimize(bb(), function(p) sum((p - c(1, 2))^2), c(0, 0),
         gr = function(p) 2 * (p - c(1, 2)))@iterations
#> [1] 2
```
