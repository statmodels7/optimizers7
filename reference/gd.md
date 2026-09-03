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
  keep_trace = FALSE
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
  Defaults to Armijo backtracking, which suffices for a method carrying
  no curvature approximation to protect.

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
repairs that by combining each direction with the previous one at no
extra cost per iteration.

Gradient descent is a baseline. On Rosenbrock from the customary start
it exhausts its budget of 500 iterations at a value of `1.9e-03` after
4905 objective evaluations, where
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md)
reaches `7.4e-10` in 25 iterations,
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md)
`5.2e-09` in 58 and
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
`3.2e-13` in 35.

## References

Cauchy's method is the oldest of them; the modern treatment, including
why its rate is linear in the condition number, is chapter 3 of Nocedal,
J. and Wright, S. J. (2006). *Numerical Optimization*, 2nd edition.
Springer, New York.

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
#>   settings  : step = 1, line_search = Armijo backtracking (c1 = 1e-04)
gd(criterion = crit_grad(1e-10), maxit = 2000)
#> <optimizer> gradient descent
#>   stop when : gradient (max-norm) < 1e-10
#>   budgets   : maxit 2000, evaluations Inf
#>   settings  : step = 1, line_search = Armijo backtracking (c1 = 1e-04)

# It solves an easy problem readily.
minimize(gd(), function(p) sum((p - c(1, 2))^2), c(0, 0),
         gr = function(p) 2 * (p - c(1, 2)))
#> <optimizer_result> gradient descent
#>   value      : 0
#>   par        : 1 2
#>   iterations : 1   evaluations: f 3, g 2
#>   elapsed    : 0 us
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)

# And is the wrong tool for a curved valley: the budget runs out first.
ros <- test_problems("rosenbrock")[[1]]
r <- minimize(gd(), ros$fn, ros$par, gr = ros$gr)
c(converged = r@converged, iterations = r@iterations,
  evaluations = r@counts[["f"]], value = r@value)
#>    converged   iterations  evaluations        value 
#> 0.000000e+00 5.000000e+02 4.905000e+03 1.915719e-03 
```
