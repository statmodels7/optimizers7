
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

# optimizers7 <img src="man/figures/logo.png" align="right" height="139" alt="" />

Almost every R package that fits a model carries its own optimiser,
written inside the function that needs it: a loop, a
`while (!converged)`, a tolerance compared against whatever quantity the
author had to hand. Nothing outside can reuse it, extend it, or ask it
for anything its author did not happen to need. The stopping rule, which
decides what the whole procedure means by *finished*, is usually a
number buried three levels down.

`{optimizers7}` writes them once, as objects. An optimiser carries its
algorithm and every setting that algorithm obeys; a stopping rule is a
separate object that can be replaced, combined, or written from scratch;
and every algorithm is written once against one objective interface.

The package is deliberately narrow. An optimiser here minimises the
function it is given and knows nothing else: not what an observation is,
not where the data came from, not what the parameters mean. Everything
that needs that knowledge belongs to the caller that has it.

It is the optimisation layer of
[statmodels7](https://statmodels7.github.io), alongside
[linkfunctions7](https://statmodels7.github.io/linkfunctions7/) and
[distributions7](https://statmodels7.github.io/distributions7/).

## Installation

``` r
# install.packages("pak")
pak::pak("statmodels7/optimizers7")
```

Or the whole toolkit at once, which also installs the four sibling
packages:

``` r
pak::pak("statmodels7/statmodels7")
```

## The shape of it

Everything minimises, always, through one generic:

``` r
f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))

minimize(bfgs(), f, par = c(-1.2, 1), gr = gr)
#> <optimizer_result> BFGS
#>   value      : 7.53493e-17
#>   par        : 1 1
#>   iterations : 32   evaluations: f 64, g 44
#>   elapsed    : 0 us
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-12 (relative))
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

A rule buried inside an algorithm fixes, at the moment the package is
written, what everyone downstream is allowed to mean by convergence.
Here it is a value:

``` r
crit_any(crit_grad(1e-10), crit_rel_obj(1e-14))
#> <criterion> gradient (max-norm) < 1e-10 or |df| < 1e-14 (relative)
```

Rules compose, and a new rule is one class and one method:

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
every point proposed is admissible by construction: there is no
rejection step, no boundary for a line search to trip over, and any
method gets bounds without knowing they exist:

``` r
y <- rnorm(200, mean = 0, sd = 2)
nll <- function(p) length(y) * log(p) + sum(y^2) / (2 * p^2)

minimize(bfgs(), nll, par = 1, lower = 0)@par
#> [1] 1.854906
sqrt(mean(y^2))
#> [1] 1.854906
```

## Starting values that need not be written out

`par` may be a starter rather than a vector. Both starters work on the
*unconstrained* scale and are mapped back through the bounds. Working on
that scale makes one constant sensible for every kind of parameter,
since zero becomes one for a variance and one half for a probability,
and it makes it impossible for a draw to land outside its box:

``` r
X <- cbind(1, matrix(rnorm(60), 30))
b <- drop(X %*% c(0.5, -1, 2)) + rnorm(30, sd = 0.2)
rss <- function(p) sum((b - X %*% p)^2)

minimize(bfgs(), rss, start_zeros())@par        # how many? worked out
#> [1]  0.5175641 -1.0649279  1.9870314
minimize(bfgs(), rss, start_runif(-2, 2))@par   # or drawn, on that scale
#> [1]  0.5175641 -1.0649279  1.9870314
```

The number of parameters is settled by whichever of three things the
caller supplies: `start_zeros(npar = 3)`, bounds with one element per
parameter, or, as above, the objective itself, probed once before the
run. `X %*% p` of the wrong length is an error rather than a number, so
a model of any kind answers the question on its own. A vectorised toy
that accepts every length cannot, and the refusal names the problem
rather than guessing.

## When the derivative does not exist

A descent method walks to the minimum of a sum of absolute deviations
and then reports that it did not converge. The report is correct,
because the subgradient it evaluates there has norm 1 and never becomes
small. `bundle()` keeps a *collection* of subgradients and tests whether
zero lies in their convex hull, which is the statement that holds at the
minimiser:

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
#> 0.06462704 0.06462703
```

## What is here

|  |  |
|----|----|
| second order | `newton()`, `bfgs()`, `lbfgs()` |
| first order | `cg()`, `bb()`, `gd()` |
| noisy objectives, and long parameter vectors | `adam()` |
| no derivative at all | `nelder_mead()`, `compass()` |
| non-smooth, with subgradients | `bundle()` |
| wrapping any of them | `multistart()`, parallel by default |
| line searches | `armijo()`, `wolfe()`, `nonmonotone()` |
| starting values | `start_zeros()`, `start_runif()` |
| stopping rules | `crit_grad()`, `crit_abs_obj()`, `crit_rel_obj()`, `crit_abs_par()`, `crit_rel_par()`, `crit_stationary()`, `crit_never()`, and `crit_any()` / `crit_all()` to combine them |

Every method reports which safeguards fired, counts its own evaluations,
and sets `converged` only when the stopping rule was satisfied, never
because the iteration budget ran out.

`multistart()` runs its starts in parallel and there is nothing to set
up: `ncores` is a number, defaulting to as many processes as there are
starts and never more than the machine can spare. Starting the workers,
loading the package on them, giving each an independent stream and
shutting them down again all happen inside and are undone before it
returns. The implementation uses forks on Linux and macOS and a socket
cluster on Windows, and `set.seed()` reproduces the run identically on
either, however many cores were used.

## Validating an optimiser

`check_optimizer()` verifies the contract: that `value` is the objective
at `par`, that a reported gradient is the gradient there, that bounds
hold *strictly*, and that a run repeats. It then reports how the method
did on the standard test problems, separately, because performance and
correctness are different questions:

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
#>     rosenbrock   gap  7.53e-17  conv       64 evals  
#>     booth        gap  4.80e-18  conv       15 evals  
#>     beale        gap  4.11e-14  conv       20 evals  
#>     powell       gap  2.86e-13  conv       60 evals  
#>     himmelblau   gap  5.54e-15  conv       18 evals  multimodal
#>     rastrigin    gap  7.96e+00  conv       15 evals  multimodal
#>     abs_sum      gap  1.26e-02  -          84 evals  non-smooth
```

## Related

- [linkfunctions7](https://statmodels7.github.io/linkfunctions7/) — link
  functions with exact derivatives to fourth order
- [distributions7](https://statmodels7.github.io/distributions7/) —
  distributions carrying exact derivatives of the log-likelihood
- [the book](https://statmodels7.github.io/book/) — the mathematics
  behind the toolkit
