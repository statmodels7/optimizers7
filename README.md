
<!-- README.md is generated from README.Rmd. Please edit that file, then
     regenerate with devtools::build_readme(). Do not use knitr::knit(): it
     processes the code but leaves this YAML header in the output as literal
     text, which GitHub and pkgdown both render verbatim. -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/statmodels7/optimizers7/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/statmodels7/optimizers7/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/statmodels7/optimizers7/graph/badge.svg)](https://app.codecov.io/gh/statmodels7/optimizers7)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# optimizers7

Almost every R package that fits a model carries its own optimiser,
written inside the function that needs it: a loop, a
`while (!converged)`, a tolerance compared against whatever quantity the
author had to hand. Nothing outside can reuse it, extend it, or ask it
for anything its author did not happen to need — and the stopping rule,
which decides what the whole thing means by *finished*, is usually a
number buried three levels down.

`{optimizers7}` writes them once, as objects. An optimiser carries its
algorithm and every setting that algorithm obeys; a stopping rule is a
separate object you can replace, combine, or write yourself; and every
algorithm is written once against one objective interface.

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

minimize(bfgs(), nll, par = 1, bounds = list(c(0, Inf)))@par
#> [1] 1.854906
sqrt(mean(y^2))
#> [1] 1.854906
```

## When the derivative does not exist

A descent method walks to the minimum of a sum of absolute deviations
and then reports that it did not converge — correctly, because the
subgradient it evaluates there has norm 1 and never becomes small.
`bundle()` keeps a *collection* of subgradients and tests whether zero
lies in their convex hull, which is the statement that actually holds:

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
| second order | `newton()`, `bfgs()`, `lbfgs()` |
| first order | `cg()`, `bb()`, `gd()` |
| noisy objectives, and long parameter vectors | `adam()` |
| no derivative at all | `nelder_mead()`, `compass()` |
| non-smooth, with subgradients | `bundle()` |
| wrapping any of them | `multistart()` |
| line searches | `armijo()`, `wolfe()` |
| stopping rules | `crit_grad()`, `crit_abs_obj()`, `crit_rel_obj()`, `crit_abs_par()`, `crit_rel_par()`, `crit_stationary()`, `crit_never()`, and `crit_any()` / `crit_all()` to combine them |

Every method reports which safeguards fired, counts its own evaluations,
and sets `converged` only when the stopping rule was satisfied — never
because the iteration budget ran out.

## Checking an optimiser of your own

`check_optimizer()` verifies the contract — that `value` is the
objective at `par`, that a reported gradient is the gradient there, that
bounds hold *strictly*, that a run repeats — and then reports how the
method did on the standard test problems, separately, because those are
different questions:

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
