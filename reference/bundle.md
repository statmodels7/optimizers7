# The Proximal Bundle Method

Proximal bundle method for convex non-smooth objectives: subgradients
collected at the visited points define a cutting-plane model, a proximal
term keeps the step near the stability centre, and convergence is tested
on the aggregate subgradient rather than on any single one, which does
not vanish at a kink.

## Usage

``` r
bundle(
  criterion = crit_stationary(1e-08),
  t0 = 1,
  t_min = 1e-10,
  t_max = 1e+10,
  m_serious = 0.1,
  bundle_size = 20,
  qp_iters = 500,
  qp_tol = 1e-12,
  maxit = 500,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 10,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule. Defaults to `crit_stationary(1e-8)`, on the
  predicted decrease.

- t0:

  Initial proximity weight, expressed as a **step length** on the
  parameter scale rather than as a bare multiplier; see Details.
  Defaults to `1`.

- t_min, t_max:

  Bounds on it, so that neither a run of null steps nor a run of serious
  ones can drive it to zero or to infinity. Default `1e-10` and `1e10`.

- m_serious:

  Fraction of the predicted decrease that a step must actually deliver
  to be accepted, in \\(0, 1)\\. Defaults to `0.1`.

- bundle_size:

  Largest number of linearisations kept before the oldest are replaced
  by their aggregate. Defaults to `20`.

- qp_iters, qp_tol:

  Effort spent on the subproblem: at most this many accelerated
  projected-gradient steps, stopping when the weights move by less than
  `qp_tol`. Defaults `500` and `1e-12`.

- maxit:

  Maximum iterations. Defaults to 500.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`: no evaluation
  budget, so the run stops on the criterion or on `maxit`. Set a finite
  value to cap the cost of a run.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 10.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

## Value

An S7 object of class `Bundle`, inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Details

### Behaviour at a kink

At the minimum of \\\lvert x \rvert\\ the evaluated subgradient is \\\pm
1\\. A descent method therefore sees a large gradient, proposes a step,
finds no acceptable one, and stops reporting failure while standing
exactly on the answer — and no tolerance can fix that, because the
quantity it is testing does not become small.

The bundle method does not test any single subgradient. It keeps a
*collection* of them, from the points it has visited, and builds the
piecewise-linear model \\\max_j \\ f_j + g_j'(x - x_j) \\\\ — the
largest of the linearisations. A kink is exactly what such a model
represents well: two linearisations meeting. At the minimum of \\\lvert
x \rvert\\ the bundle holds subgradients near \\+1\\ and near \\-1\\,
and \\0\\ is a convex combination of them, which is precisely the
statement that \\0\\ lies in the subdifferential.

That is what is tested. The step solves \$\$\min_d\\ \max_j \\
-\alpha_j + g_j'd \\ + \frac{1}{2t}\lVert d \rVert^2,\$\$ where
\\\alpha_j \ge 0\\ measures how badly linearisation \\j\\ misses the
current point and the quadratic term keeps the step inside the region
where the model is believed. Its optimal value \\v \le 0\\ is the
*predicted decrease*, which is what the acceptance test below uses.

What
[`crit_stationary`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
watches is a related but different quantity, the *optimality estimate*
\\\lVert p \rVert^2 + \alpha\\, where \\p\\ is the aggregate
subgradient. Both vanish together at a solution, and the distinction is
not pedantry: \\-v\\ carries a factor of \\t\\, and \\t\\ is halved at
every null step, so it can be pushed below any tolerance by the trust
parameter shrinking rather than by the point becoming stationary. A run
doing that reports success while standing somewhere its own model still
calls steep. Dropping the factor leaves the model's own claim, which no
amount of shrinking can flatter.

### Serious steps and null steps

A trial point is accepted — a *serious step* — only if it delivers at
least `m_serious` of the decrease the model promised. That is the same
bargain an Armijo condition strikes, for the same reason: accepting any
decrease at all lets a sequence of ever tinier improvements masquerade
as progress.

When it is refused, the iteration is not wasted. The trial point
contributes its subgradient to the bundle, so the model is strictly
better next time; this is a *null step*, and a run that reports many of
them is refining its picture of a kink, not failing. Both counts appear
in the result's message.

`t` is halved after a null step and doubled after a serious one, within
`t_min` and `t_max`. Kiwiel's rule chooses the factor from a curvature
estimate and is better; this is the safeguarded version, which is cruder
and bounded.

### The trust parameter t0

The step is \\d = -t\\p\\, so a bare `t` would make the first step as
long as the gradient happens to be big. On Rosenbrock from its usual
start that is 233 units, landing where the objective is \\10^{11}\\ and
its gradient \\10^{15}\\; the null step that follows halves `t` while
the gradient has squared, so `t` can never catch up and the run spends
its entire budget on rejected steps before returning the point it began
at.

`t0` is therefore divided by the size of the gradient at the starting
point, which makes the first step of length `t0` in parameter space. It
is the same normalisation a line search performs when it starts at 1
along a unit direction, and `t_min` and `t_max` move with it since they
bound the same quantity. Raise `t0` for a problem whose optimum is far
away, lower it for one where the model is trustworthy only nearby.

Should a subgradient overflow anyway — possible for an objective that
grows fast enough — the run stops and says so, rather than spending its
budget on a subproblem whose matrix contains an infinity.

### Bounded memory

Left alone the bundle grows without limit. When it reaches `bundle_size`
the oldest linearisations are replaced by the *aggregate* — the single
affine function the subproblem's solution defines — rather than
discarded. Discarding loses what they knew and can stall the method; the
aggregate keeps a summary of all of it in one element, which is what
makes a bounded bundle safe.

### Convexity requirement

The theory is for convex \\f\\, where the linearisation errors
\\\alpha_j\\ are non-negative automatically. On a non-convex objective
they can come out negative and are clipped at zero. The method then
still runs and usually works, but that clip is exactly where the
guarantee stops: a negative error is the model reporting that \\f\\ lies
below its own linearisation, and setting it to zero suppresses the
information rather than using it.

### Subgradients

Supply `gr`. Without one the package differences the objective, and
although a difference is a perfectly good gradient wherever \\f\\ is
differentiable — which is almost everywhere, so a run mostly gets away
with it — a difference taken *across* a kink is not a subgradient of
anything and the model will be built from a number that belongs to no
linearisation.

## References

Kiwiel, K. C. (1990). Proximity control in bundle methods for convex
nondifferentiable minimization. *Mathematical Programming* **46**,
105–122.

Mäkelä, M. M. (2002). Survey of bundle methods for nonsmooth
optimization. *Optimization Methods and Software* **17**, 1–29.

## See also

[`nelder_mead`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
[`compass`](https://statmodels7.github.io/optimizers7/reference/compass.md),
[`crit_stationary`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)

## Examples

``` r
bundle()
#> <optimizer> proximal bundle
#>   stop when : stationarity < 1e-08
#>   budgets   : maxit 500, evaluations Inf
#>   settings  : t0 = 1, t_min = 1e-10, t_max = 1e+10, m_serious = 0.1, bundle_size = 20, qp_iters = 500, qp_tol = 1e-12

# the median, as the minimiser of a sum of absolute deviations: the objective
# has a kink at every observation, and one of them is the answer
set.seed(1)
y <- rnorm(101)
f <- function(p) sum(abs(y - p))
g <- function(p) -sum(sign(y - p))
r <- minimize(bundle(), f, par = 0, gr = g)
c(bundle = r@par, median = median(y))
#>     bundle     median 
#> 0.07456498 0.07456498 

# a kink running diagonally, where a coordinate-wise search stalls
minimize(bundle(), function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2),
         c(1, 0.5), gr = function(p) c(sign(p[1] + p[2]), sign(p[1] + p[2])) +
                         0.2 * p)@value
#> [1] 8.207375e-09
```
