# optimizers7

Almost every R package that fits a model carries its own optimizer,
written inside the function that needs it: a loop, a
`while (!converged)`, a tolerance compared with a quantity chosen by the
author. Such an optimizer cannot be reused or extended outside the
package, and it provides only the quantities that its author needed. The
stopping rule, which defines what the procedure means by *finished*, is
usually a constant deep inside the code.

In [optimizers7](https://statmodels7.github.io/optimizers7/), each
optimizer is written once, as an object. An optimizer carries its
algorithm and every setting of that algorithm; a stopping rule is a
separate object that can be replaced, combined, or written by the user;
and every algorithm is written once against one interface for the
objective.

The package minimizes the function that it is given, and it does not
represent observations, the origin of the data or the meaning of the
parameters. Those belong to the calling package.

It is the optimization layer of
[statmodels7](https://statmodels7.github.io), alongside
[linkfunctions7](https://statmodels7.github.io/linkfunctions7/) and
[distributions7](https://statmodels7.github.io/distributions7/).

## Installation

``` r

# install.packages("pak")
pak::pak("statmodels7/optimizers7")
```

To install the whole toolkit, including the other packages of the
toolkit, use:

``` r

pak::pak("statmodels7/statmodels7")
```

## Minimizing a function

Every algorithm minimizes, through one generic:

``` r

f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))

minimize(bfgs(), f, par = c(-1.2, 1), gr = gr)
#> <optimizer_result> BFGS
#>   value      : 3.21085e-13
#>   par        : 1 1
#>   iterations : 35   evaluations: f 49, g 39
#>   elapsed    : 0 us
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
```

Changing the algorithm changes one word, because every optimizer carries
the same shared settings:

``` r

minimize(lbfgs(), f, c(-1.2, 1), gr = gr)@par
#> [1] 1 1
minimize(newton(), f, c(-1.2, 1), gr = gr)@par
#> [1] 1 1
minimize(nelder_mead(), f, c(-1.2, 1))@par
#> [1] 1 1
```

## Stopping rules

A rule written inside an algorithm fixes the meaning of convergence when
the package is written. Here a rule is an object:

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

A rule that an optimizer cannot evaluate is rejected when the run
starts, with an error that names the missing quantity, instead of being
accepted and never firing:

``` r

minimize(nelder_mead(criterion = crit_grad()), f, c(-1.2, 1))
#> Error:
#> ! The stopping rule needs gradient, which nelder-mead does not provide.
#>   Choose a criterion this optimizer can evaluate, or a method that provides it.
```

## Bounds by reparametrization

Each bounded coordinate is reparametrized onto the whole real line, so
every proposed point is admissible by construction: there is no
rejection step and no boundary for a line search to cross, and a method
needs no code of its own for bounds. The exception is
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md),
which takes its constraint through its proximal operator and rejects
`lower` and `upper`:

``` r

y <- rnorm(200, mean = 0, sd = 2)
nll <- function(p) length(y) * log(p) + sum(y^2) / (2 * p^2)

minimize(bfgs(), nll, par = 1, lower = 0)@par
#> [1] 1.854906
sqrt(mean(y^2))
#> [1] 1.854906
```

## Starting values that need not be written out

`par` may be a starter instead of a vector. Both starters work on the
*unconstrained* scale and are mapped back through the bounds. On that
scale one constant suits every kind of parameter, since zero becomes one
for a variance and one half for a probability, and every draw lies
inside its box:

``` r

X <- cbind(1, matrix(rnorm(60), 30))
b <- drop(X %*% c(0.5, -1, 2)) + rnorm(30, sd = 0.2)
rss <- function(p) sum((b - X %*% p)^2)

minimize(bfgs(), rss, start_zeros())@par        # how many? worked out
#> [1]  0.5175641 -1.0649279  1.9870314
minimize(bfgs(), rss, start_runif(-2, 2))@par   # or drawn, on that scale
#> [1]  0.5175641 -1.0649279  1.9870314
```

The number of parameters is given by `start_zeros(npar = 3)`, by bounds
with one element per parameter, or, as above, by the objective itself,
probed once before the run. `X %*% p` with a parameter of the wrong
length is an error and not a number, so an objective written as a matrix
product determines the length. A vectorized expression that accepts
several lengths does not, and the probe then raises an error that names
the lengths accepted.

## Objectives without a derivative at the minimum

A descent method reaches the minimum of a sum of absolute deviations and
then reports that it did not converge. The report is correct, because
the subgradient that it evaluates near the minimum has norm 1 and does
not become small.
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
keeps a *collection* of subgradients and tests whether zero lies in
their convex hull, which is the condition that holds at the minimizer:

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
#> 0.06462703 0.06462703
```

## Components of the package

|  |  |
|----|----|
| second order | [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md), [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md), [`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md) |
| third order | [`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md) |
| first order | [`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md), [`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md), [`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md) |
| noisy objectives, and long parameter vectors | [`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md) |
| no derivative | [`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md), [`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md) |
| one parameter | [`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md) |
| global search | [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) |
| non-smooth, with subgradients | [`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md) |
| non-smooth, with a proximal operator | [`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md) |
| wrappers | [`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md), parallel by default, and [`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md) |
| line searches | [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md), [`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md), [`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md) |
| starting values | [`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md), [`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md) |
| stopping rules | [`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md), [`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md), [`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md), [`crit_abs_par()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md), [`crit_rel_par()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_par.md), [`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md), [`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md), and [`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md) / [`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md) to combine them |

Every method counts its own evaluations, records in its trace the
safeguards that fired, and sets `converged` only when the run ended on
its stopping rule (or on the resolution of the objective), never because
a budget ran out.

[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
runs its starts in parallel without any setup: `cores` is a number, by
default as many processes as there are starts and at most the number of
cores minus two. The call itself starts the workers, loads the package
on them, gives each of them an independent random stream and shuts them
down. The implementation uses forks on Linux and macOS and a socket
cluster on Windows. For a deterministic optimizer such as
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
[`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces the run
identically however many cores were used; for a stochastic one such as
[`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md) the
result depends on `cores`.

## Validating an optimizer

[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
runs twelve checks on what an optimizer reports: that `value` is the
objective at `par`, that a reported gradient is the gradient there, that
bounds hold *strictly*, that a run repeats, and others. It then reports,
separately, the gaps that the method reached on the standard test
problems; these gaps measure performance and do not enter the verdicts:

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
#>   [10] an unevaluable rule is rejected:  [PASSED]
#>   [11] a bad starting point is an error: [PASSED]
#>   [12] it minimizes a quadratic:         [PASSED]
#> 
#>   All checks passed.
#> 
#>   battery (gap from the known minimum; information, not a verdict)
#>     sphere       gap  4.62e-33  conv        3 evals  
#>     rosenbrock   gap  3.21e-13  conv       49 evals  
#>     booth        gap  1.73e-17  conv        8 evals  
#>     beale        gap  3.18e-15  conv       17 evals  
#>     powell       gap  4.28e-11  conv       39 evals  
#>     himmelblau   gap  6.27e-14  conv       17 evals  multimodal
#>     rastrigin    gap  0.00e+00  conv       11 evals  multimodal
#>     abs_sum      gap  1.26e-02  -          84 evals  non-smooth
```

## Related

- [linkfunctions7](https://statmodels7.github.io/linkfunctions7/): link
  functions with exact derivatives to fifth order
- [distributions7](https://statmodels7.github.io/distributions7/):
  distributions with exact derivatives of the log-likelihood
- [the book](https://statmodels7.github.io/book/): the mathematics of
  the toolkit
