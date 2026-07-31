# Gradient Descent

The simplest method there is: step along the negative gradient, and let
the line search decide how far.

## Usage

``` r
gd(
  criterion = crit_any(crit_grad(1e-08), crit_rel_obj(1e-12)),
  step = 1,
  line_search = armijo(),
  maxit = 500,
  max_eval = 10000,
  verbose = FALSE,
  refresh = 10,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule; see
  [`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md).

- step:

  Initial step length offered to the line search each iteration.
  Defaults to `1`.

- line_search:

  How far to go along the direction; see
  [`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
  and
  [`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md).
  Defaults to Armijo backtracking, which is all a method with nothing to
  protect needs.

- maxit:

  Maximum iterations. Defaults to 500.

- max_eval:

  Maximum objective evaluations. Defaults to 10000.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 10.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

## Value

An S7 object of class `GradientDescent`, inheriting from
[`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Details

The direction \\-g\\ minimises \\g^\top d\\ over directions of a given
length, so it is the steepest one available — and that is also the whole
criticism of it, because steepest is a statement about an arbitrary
choice of length. Consecutive directions under an exact line search are
orthogonal, so on a narrow valley the iterates zigzag across the floor
and progress along it only slowly.
[`cg`](https://statmodels7.github.io/optimizers7/reference/cg.md)
repairs exactly this by bending each direction back towards the last, at
no extra cost per iteration.

Reach for it when you want the simplest thing that works, or as the
baseline another method has to beat. On anything smooth
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/Bfgs.md)
will beat it by orders of magnitude.

## See also

[`cg`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`bb`](https://statmodels7.github.io/optimizers7/reference/bb.md),
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/Bfgs.md)

## Examples

``` r
gd()
#> <optimizer> gradient descent
#>   stop when : gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
#>   budgets   : maxit 500, evaluations 10000
#>   settings  : step = 1, line_search = Armijo backtracking (c1 = 1e-04)
gd(criterion = crit_grad(1e-10), maxit = 2000)
#> <optimizer> gradient descent
#>   stop when : gradient (max-norm) < 1e-10
#>   budgets   : maxit 2000, evaluations 10000
#>   settings  : step = 1, line_search = Armijo backtracking (c1 = 1e-04)

minimize(gd(), function(p) sum((p - c(1, 2))^2), c(0, 0),
         gr = function(p) 2 * (p - c(1, 2)))
#> <optimizer_result> gradient descent
#>   value      : 0
#>   par        : 1 2
#>   iterations : 1   evaluations: f 3, g 2
#>   converged  : yes (gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative))
```
