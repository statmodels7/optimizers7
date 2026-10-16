# Line Search Satisfying the Strong Wolfe Conditions

Searches for a step that satisfies both Wolfe conditions, sufficient
decrease and curvature: \$\$f(x + \alpha d) \leq f(x) + c_1 \alpha\\
g^\top d, \qquad \lvert g(x + \alpha d)^\top d \rvert \leq c_2 \lvert
g^\top d \rvert .\$\$ It evaluates the gradient at trial points as well
as the objective, which is the cost of the second condition, the one
that
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
need.

## Usage

``` r
wolfe(c1 = 1e-04, c2 = 0.9, max_step = 30, resolution = 0)
```

## Arguments

- c1:

  Sufficient-decrease constant, strictly inside \\(0, 1)\\. Defaults to
  `1e-4`.

- c2:

  Curvature constant, with \\c_1 \< c_2 \< 1\\. Defaults to `0.9`, the
  usual choice for a quasi-Newton method;
  [`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md)
  uses `0.1`, a more exact line search. A `c2` at or below `c1` is
  rejected, because a step that satisfies both conditions is guaranteed
  to exist only when \\c_1 \< c_2 \< 1\\.

- max_step:

  Maximum trial steps in each of the bracketing and zoom phases, a
  positive whole number. Defaults to 30.

- resolution:

  The smallest difference in the objective that the objective can
  resolve, in the units of the objective, or a function of no arguments
  that returns it when it changes during the run. Defaults to `0`, which
  disables the test; see
  [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
  for the purpose of the test and for the reason why it is applied at
  the full step, before the search begins.

## Value

An S7 object of class
[WolfeSearch](https://statmodels7.github.io/optimizers7/reference/WolfeSearch-class.md),
inheriting from
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md),
to be passed as an optimizer's `line_search`.

## The curvature condition

[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
does not enforce it, and a quasi-Newton method relies on it. BFGS builds
its approximation from the secant pair \\(s, y)\\ with \\s = \alpha d\\
and \\y = g\_{new} - g\_{old}\\, and a step so short that the gradient
has barely changed gives a pair with no usable curvature: the update is
then skipped, or it degrades the approximation. Requiring the
directional derivative to have shrunk by the factor \\c_2\\ rules this
out. On the problems of
[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md),
BFGS under
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
skips an update on a few problems only; see
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md).

## The implementation

The search brackets a step, doubling the trial step until the bracket is
found, and then bisects inside the bracket instead of interpolating a
polynomial. Bisection costs a few more evaluations and does not depend
on the shape of the function inside the bracket. Each phase is limited
to `max_step` trials, and a search that exhausts them returns without a
step. When the bracket collapses before the curvature condition holds,
the end of the bracket that satisfies sufficient decrease is returned
and the trace records `curvature condition not met`.

The gradient evaluations at the trial points make this search more
costly per iteration than
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md);
the total cost of a run depends on the problem.

## References

Wolfe, P. (1969). Convergence conditions for ascent methods. *SIAM
Review* **11**, 226–235.

Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*, 2nd
edition. Springer, New York.

## See also

[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
for the cheap search,
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
for the relaxed reference value,
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
for the method that needs this one.

## Examples

``` r
wolfe()
#> <line_search> strong Wolfe (c1 = 1e-04, c2 = 0.9)
wolfe(c2 = 0.1)
#> <line_search> strong Wolfe (c1 = 1e-04, c2 = 0.1)

# Gradient descent with this search comes close to the minimum (1, 1) in
# 2000 iterations.
minimize(gd(line_search = wolfe(), maxit = 2000),
         function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2,
         c(-1.2, 1))@par
#> [1] 1.001140 1.002291

# The smaller constant matters for conjugate gradients, which has no
# curvature approximation to compensate for an inexact step.
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))
c(tight = minimize(cg(line_search = wolfe(c2 = 0.1)), f, c(-1.2, 1),
                   gr = g)@iterations,
  loose = minimize(cg(line_search = wolfe(c2 = 0.9)), f, c(-1.2, 1),
                   gr = g)@iterations)
#> tight loose 
#>    25    28 

# c2 must exceed c1.
try(wolfe(c1 = 0.5, c2 = 0.4))
#> Error : 'c2' must be greater than 'c1'; only then is a step that satisfies both conditions guaranteed to exist.
```
