# Line Search Satisfying the Strong Wolfe Conditions

Finds a step satisfying both Wolfe conditions, sufficient decrease and
curvature together: \$\$f(x + \alpha d) \leq f(x) + c_1 \alpha\\ g^\top
d, \qquad \lvert g(x + \alpha d)^\top d \rvert \leq c_2 \lvert g^\top d
\rvert .\$\$ Evaluates the gradient at trial points as well as the
objective, which is what the second condition costs and what
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
  defaults to `0.1`, wanting a more exact line search. A `c2` at or
  below `c1` is refused, the two conditions being unsatisfiable together
  then.

- max_step:

  Maximum trial steps in **each** of the bracketing and zoom phases, a
  positive whole number. Defaults to 30.

- resolution:

  The smallest difference in the objective that means anything, in the
  objective's own units, or a function of no arguments returning it
  where it moves as the run goes. Defaults to `0`, which does not ask
  the question; see
  [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
  for what it is for and why the question is put at the full step,
  before the search begins.

## Value

An S7 object of class
[WolfeSearch](https://statmodels7.github.io/optimizers7/reference/WolfeSearch-class.md),
inheriting from
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md),
to be passed as an optimizer's `line_search`.

## What the curvature condition is for

[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
cannot provide it, and a quasi-Newton method needs it to work at all.
BFGS builds its approximation from the secant pair \\(s, y)\\ with \\s =
\alpha d\\ and \\y = g\_{new} - g\_{old}\\, and a step so short that the
gradient has barely moved gives a pair carrying no curvature: the update
is then either skipped or it corrupts the matrix. Requiring the gradient
along the direction to have shrunk by a factor \\c_2\\ is exactly the
guarantee that this does not happen.

In practice the guarantee is worth less than it sounds, and the honest
measurement belongs here: over the eight problems of
[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md),
BFGS under Armijo skips an update on three of them, once each. See
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
for the numbers.

## The implementation

Bracketing and zoom, with **bisection** inside the zoom instead of
polynomial interpolation. That costs a few more evaluations and cannot
be defeated by an awkwardly shaped interval.

Gradient evaluations at trial points are the price, so this is the
dearer choice per iteration and usually the cheaper one per problem.

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

# Even gradient descent gets to the answer with it, given the budget.
minimize(gd(line_search = wolfe(), maxit = 2000),
         function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2,
         c(-1.2, 1))@par
#> [1] 1.001140 1.002291

# The tighter constant matters to conjugate gradients, which has no
# curvature approximation to repair a loose step with.
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))
c(tight = minimize(cg(line_search = wolfe(c2 = 0.1)), f, c(-1.2, 1),
                   gr = g)@iterations,
  loose = minimize(cg(line_search = wolfe(c2 = 0.9)), f, c(-1.2, 1),
                   gr = g)@iterations)
#> tight loose 
#>    25    28 

# The two conditions have to be satisfiable together.
try(wolfe(c1 = 0.5, c2 = 0.4))
#> Error : 'c2' must be greater than 'c1'; the conditions are unsatisfiable otherwise.
```
