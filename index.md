# optimizers7

Almost every R package that fits a model carries its own optimiser,
written inside the function that needs it: a loop, a
`while (!converged)`, a tolerance compared against whatever quantity the
author had to hand. Nothing outside can reuse it, extend it, or ask it
for anything its author did not happen to need — and the stopping rule,
which decides what the whole thing means by *finished*, is usually a
number buried three levels down.

[optimizers7](https://statmodels7.github.io/optimizers7/) writes them
once, as objects. An optimiser carries its algorithm and every setting
that algorithm obeys; a stopping rule is a separate object you can
replace, combine, or write yourself; and every algorithm is written once
against one objective interface.

It is deliberately narrow. An optimiser here minimises the function it
is given and knows nothing else — not what an observation is, not where
the data came from, not what the parameters mean. Everything that needs
that knowledge belongs to the caller that has it.

It is the optimisation layer of
[statmodels7](https://statmodels7.github.io), alongside
[linkfunctions7](https://statmodels7.github.io/linkfunctions7) and
[distributions7](https://statmodels7.github.io/distributions7).

## Installation

``` r

# install.packages("pak")
pak::pak("statmodels7/optimizers7")
```

## The shape of it

Everything minimises, always, through one generic:

``` r

f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))

minimize(bfgs(), f, par = c(-1.2, 1), gr = gr)
#> <optimizer_result> BFGS
#>   value      : 2.04447e-20
#>   par        : 1 1
#>   iterations : 33   evaluations: f 65, g 45
#>   converged  : yes (gradient (max-norm) < 1e-08 or |df| < 1e-12 (relative))
```

Swapping the algorithm changes one word, because every optimiser carries
the same settings:

``` r

minimize(lbfgs(), f, c(-1.2, 1), gr = gr)@par
#> [1] 1 1
minimize(newton(), f, c(-1.2, 1), gr = gr)@par
#> [1] 1 1
minimize(nelder_mead(), f, c(-1.2, 1))@par
#> [1] 1 1
```

## The stopping rule is an object

This is the part that is not a convenience. A rule buried inside an
algorithm fixes, at the moment the package is written, what everyone
downstream is allowed to mean by convergence. Here it is a value:

``` r

crit_any(crit_grad(1e-10), crit_rel_obj(1e-14))
#> <criterion> gradient (max-norm) < 1e-10 or |df| < 1e-14 (relative)
```

Rules compose, and a rule the authors never thought of is a first-class
citizen — one class, one method:

``` r

CritTiny <- S7::new_class("CritTiny", parent = criterion,
                          properties = list(tol = S7::class_numeric))
S7::method(crit_met, CritTiny) <- function(criterion, state)
  state$f_new < criterion@tol

minimize(bfgs(criterion = CritTiny(label = "f < 1e-8", tol = 1e-8)),
         f, c(-1.2, 1), gr = gr)@criterion_met
#> [1] "f < 1e-8"
```

A rule an optimiser cannot evaluate is **refused**, by name, rather than
accepted and left never to fire:

``` r

minimize(nelder_mead(criterion = crit_grad()), f, c(-1.2, 1))
#> Error:
#> ! The stopping rule needs gradient, which nelder-mead does not provide.
#>   Choose a criterion this optimiser can evaluate, or a method that provides it.
```

## Bounds are removed, not enforced

Each bounded coordinate is reparametrised onto the whole real line, so
every point proposed is admissible by construction — no rejection step,
no boundary for a line search to trip over, and any method gets bounds
without knowing they exist:

``` r

y <- rnorm(200, mean = 0, sd = 2)
nll <- function(p) length(y) * log(p) + sum(y^2) / (2 * p^2)

minimize(bfgs(), nll, par = 1, lower = 0)@par
#> [1] 1.854906
sqrt(mean(y^2))
#> [1] 1.854906
```

## When the derivative does not exist

A descent method walks to the minimum of a sum of absolute deviations
and then reports that it did not converge — correctly, because the
subgradient it evaluates there has norm 1 and never becomes small.
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/Bundle.md)
keeps a *collection* of subgradients and tests whether zero lies in
their convex hull, which is the statement that actually holds:

``` r

z <- rnorm(101)
sad <- function(p) sum(abs(z - p))
sub <- function(p) -sum(sign(z - p))

c(bfgs   = minimize(bfgs(), sad, par = 0, gr = sub)@converged,
  bundle = minimize(bundle(), sad, par = 0, gr = sub)@converged)
#>   bfgs bundle 
#>  FALSE   TRUE

c(bundle = minimize(bundle(), sad, par = 0, gr = sub)@par, median = median(z))
#>     bundle     median 
#> 0.00213186 0.00213186
```

## What is here

|  |  |
|----|----|
| second order | [`newton()`](https://statmodels7.github.io/optimizers7/reference/Newton.md), [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/Bfgs.md), [`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/Lbfgs.md) |
| first order | [`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md), [`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md), [`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) |
| noisy objectives, and long parameter vectors | [`adam()`](https://statmodels7.github.io/optimizers7/reference/Adam.md) |
| no derivative at all | [`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md), [`compass()`](https://statmodels7.github.io/optimizers7/reference/Compass.md) |
| non-smooth, with subgradients | [`bundle()`](https://statmodels7.github.io/optimizers7/reference/Bundle.md) |
| wrapping any of them | [`multistart()`](https://statmodels7.github.io/optimizers7/reference/MultiStart.md) |
| line searches | [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md), [`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md) |
| stopping rules | [`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md), [`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md), [`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md), [`crit_abs_par()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md), [`crit_rel_par()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_par.md), [`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md), [`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md), and [`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md) / [`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md) to combine them |

Every method reports which safeguards fired, counts its own evaluations,
and sets `converged` only when the stopping rule was satisfied — never
because the iteration budget ran out.

## Checking an optimiser of your own

[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
verifies the contract — that `value` is the objective at `par`, that a
reported gradient is the gradient there, that bounds hold *strictly*,
that a run repeats — and then reports how the method did on the standard
test problems, separately, because those are different questions:

``` r

check_optimizer(bfgs())
#> Checking optimizer: BFGS
#>   [ 1] value agrees with par:            [PASSED]
#>   [ 2] gradient agrees with par:         [PASSED]
#>   [ 3] convergence is not assumed:       [PASSED]
#>   [ 4] budgets are respected:            [PASSED]
#>   [ 5] evaluations are counted:          [PASSED]
#>   [ 6] trace is well formed:             [PASSED]
#>   [ 7] bounds are respected strictly:    [PASSED]
#>   [ 8] the run repeats:                  [PASSED]
#>   [ 9] maximize mirrors minimize:        [PASSED]
#>   [10] an unevaluable rule is refused:   [PASSED]
#>   [11] a bad starting point is an error: [PASSED]
#>   [12] it minimises a quadratic:         [PASSED]
#> 
#>   All checks passed.
#> 
#>   battery (gap from the known minimum; information, not a verdict)
#>     sphere       gap  0.00e+00  conv        3 evals  
#>     rosenbrock   gap  2.04e-20  conv       65 evals  
#>     booth        gap  4.80e-18  conv       15 evals  
#>     beale        gap  2.40e-17  conv       21 evals  
#>     powell       gap  3.65e-16  conv       75 evals  
#>     himmelblau   gap  3.04e-19  conv       19 evals  multimodal
#>     rastrigin    gap  7.96e+00  conv       15 evals  multimodal
#>     abs_sum      gap  1.26e-02  -          84 evals  non-smooth
```

## Related

- [linkfunctions7](https://statmodels7.github.io/linkfunctions7) — link
  functions with exact derivatives to fourth order
- [distributions7](https://statmodels7.github.io/distributions7) —
  distributions carrying exact derivatives of the log-likelihood
- [the book](https://statmodels7.github.io/book/) — the mathematics,
  derived
