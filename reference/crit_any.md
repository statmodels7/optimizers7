# Stop When Any of Several Rules Fires

Combines criteria disjunctively. This is the usual arrangement: a run
should end as soon as any reasonable rule is satisfied.

## Usage

``` r
crit_any(...)
```

## Arguments

- ...:

  [`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  objects.

## Value

A
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
object, so combinations nest.

## Details

**The default rule of the gradient methods.**
[`gd`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`cg`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`bb`](https://statmodels7.github.io/optimizers7/reference/bb.md),
[`newton`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
and
[`lbfgs`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
default to `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`: the
point is stationary, or the objective has stopped moving, or the
parameters have. Since a disjunction can only get weaker as terms are
added, a run that ends under this rule would have ended under a gradient
rule alone at best later and never earlier.

What that buys and what it costs was measured over the package's own
[`test_problems`](https://statmodels7.github.io/optimizers7/reference/test_problems.md),
six methods on eight problems. Against the gradient rule alone it
converges on 44 of the 48 runs rather than 41, and costs 19370 objective
evaluations rather than 22299. The three it gains are `cg` and `bb` on
the non-smooth `abs_sum` and `gd` on Beale, and none is lost.

The cost is that a run stops sooner, so the point it reports is further
from the solution. Measured, 10 of the 48 end at a gradient more than a
hundred times larger, and the worst are runs that were reaching absurd
precision anyway: `bfgs` on Rosenbrock ends at 8.1e-06 rather than
4.4e-10, with the objective 3.2e-13 above its minimum rather than
1.2e-21. On one run the difference is real rather than cosmetic – `cg`
on `abs_sum`, where the objective ends 4.5e-02 above the minimum rather
than 1.9e-03, and the flag reads `TRUE` where it used to read `FALSE`.
That is a smooth method on a non-smooth problem, where the objective
stalls far from the solution and a rule that reads a stall cannot tell
the two apart. A caller who needs stationarity asks for it:
`criterion = crit_grad()`.

[`crit_rel_obj`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
was in that default until 0.6.0 and is not any more, because it never
fired: measured over the same 48 runs, the rule with it and the rule
without it agree on every count, every evaluation and every reported
point. It remains available and useful where an objective's scale is not
known in advance.

## See also

[`crit_all`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)

## Examples

``` r
crit_any(crit_grad(1e-8), crit_rel_obj(1e-12))
#> <criterion> gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative)
```
