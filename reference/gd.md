# Gradient Descent

Steepest descent: the search direction is the negative gradient
\\-g(x)\\, and a line search chooses the step length.

## Usage

``` r
gd(
  criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
  step = 1,
  line_search = armijo(),
  maxit = 500,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 10,
  keep_trace = FALSE,
  stop_decrement = 0.001
)
```

## Arguments

- criterion:

  The stopping rule; see
  [`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md).

- step:

  Initial step length offered to the line search each iteration.
  Defaults to `1`.

- line_search:

  How far to go along the direction; see
  [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
  and
  [`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md).
  Defaults to Armijo backtracking, which is sufficient for a method that
  keeps no curvature approximation.

- maxit:

  Maximum iterations. Defaults to 500.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`, which sets no
  evaluation budget, so the run stops on the criterion or on `maxit`.
  Set a finite value to cap the cost of a run.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 10.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

- stop_decrement:

  The largest decrease \\\tfrac12\lvert g^\top s\rvert\\, with \\s\\ the
  last step and \\g\\ the gradient at the new point, at which a stopping
  rule on the change in the parameters or in the objective, such as
  [`crit_abs_par()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md)
  or
  [`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md),
  may end the run as converged. A single positive number, in the units
  of the objective, or `Inf` to accept every such stop. Defaults to
  `1e-3`. See below.

## Value

An S7 object of class
[GradientDescent](https://statmodels7.github.io/optimizers7/reference/GradientDescent-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Details

The direction \\-g\\ minimizes \\g^\top d\\ over directions of a given
Euclidean length. Under an exact line search consecutive directions are
orthogonal, so on an ill-conditioned objective the iterates zigzag and
converge slowly.
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md)
corrects this by combining each direction with the previous one, at no
extra cost per iteration.

## A stop on a change

A stopping rule on the change in the parameters or in the objective
states that the iterate no longer moves, and near a minimum that is the
same as being there. Away from one it is not: a step can be short
because the objective is steep, or because the direction came from a
poor model of the surface, and the decrease \\\tfrac12\lvert g^\top
s\rvert\\ that the step still predicts at the new point is then above
what a minimum leaves. A stop on the change is therefore not accepted
where that decrease exceeds `stop_decrement`. The run ends as not
converged, with `criterion_met` reporting `"stopped without converging"`
and the message giving the decrease, and the trace records the safeguard
`change stop refused`. A rule on the gradient, such as
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md),
is not affected. The check uses quantities that the loop already holds.
`stop_decrement = Inf` accepts every stop on a change.

Gradient descent serves as a baseline for the other methods. On the
Rosenbrock function from the usual start it does not converge within its
default budget of 500 iterations, while
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) and
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
converge in far fewer (the last example).

## References

Cauchy's method is the oldest of the methods in the package; chapter 3
of the following book treats it, including the dependence of its linear
rate on the condition number: Nocedal, J. and Wright, S. J. (2006).
*Numerical Optimization*, 2nd edition. Springer, New York.

## See also

[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md) for
the same storage and a better direction,
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) for
a scalar curvature estimate,
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
for a matrix one.

## Examples

``` r
gd()
#> <optimizer> gradient descent
#>   stop when : gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08
#>   budgets   : maxit 500, evaluations Inf
#>   settings  : step = 1, line_search = Armijo backtracking (c1 = 1e-04), stop_decrement = 0.001
gd(criterion = crit_grad(1e-10), maxit = 2000)
#> <optimizer> gradient descent
#>   stop when : gradient (max-norm) < 1e-10
#>   budgets   : maxit 2000, evaluations Inf
#>   settings  : step = 1, line_search = Armijo backtracking (c1 = 1e-04), stop_decrement = 0.001

# A well-conditioned problem is solved quickly.
minimize(gd(), function(p) sum((p - c(1, 2))^2), c(0, 0),
         gr = function(p) 2 * (p - c(1, 2)))
#> <optimizer_result> gradient descent
#>   value      : 0
#>   par        : 1 2
#>   iterations : 1   evaluations: f 3, g 2
#>   elapsed    : 1 ms
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)

# In a curved valley the budget runs out before convergence.
ros <- test_problems("rosenbrock")[[1]]
r <- minimize(gd(), ros$fn, ros$par, gr = ros$gr)
c(converged = r@converged, iterations = r@iterations,
  evaluations = r@counts[["f"]], value = r@value)
#>    converged   iterations  evaluations        value 
#> 0.000000e+00 5.000000e+02 4.905000e+03 1.915719e-03 
```
