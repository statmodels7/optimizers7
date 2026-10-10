# Stop When Any of Several Rules Fires

Combines criteria disjunctively: the run ends as soon as any of the
rules is satisfied.

## Usage

``` r
crit_any(...)
```

## Arguments

- ...:

  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  objects.

## Value

A
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
object, so combinations nest.

## Details

**The default rule of the gradient methods.**
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md),
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
default to `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`: the
point is stationary, or the objective has stopped moving, or the
parameters have. A disjunction can only become weaker as terms are
added, so a run that ends under this rule would have ended under the
gradient rule alone at the same iteration or later, never earlier.

On the package's
[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md),
with the six methods other than
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md),
this rule converges on every run on which the gradient rule alone
converges, and on three more (`cg` and `bb` on `abs_sum`, `gd` on
`beale`), and it needs fewer evaluations of the objective and the
gradient in total. The cost is that a run stops sooner, at a point
further from the solution, usually by a margin far below the statistical
precision of an estimate. On a non-smooth problem a smooth method can
stall far from the solution, and a rule that reads a stall reports it as
convergence (`cg` on `abs_sum` is such a case). A caller who needs a
stationary point uses `criterion = crit_grad()`.

[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
was part of this default until 0.6.0. It was removed because on the same
problems it never fired before the other terms.

## See also

[`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)

## Examples

``` r
crit_any(crit_grad(1e-8), crit_rel_obj(1e-12))
#> <criterion> gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
```
