# The Standard Test Problems

Returns the eight functions optimization papers are argued over: a
quadratic, a curved valley, a few awkward polynomials, two with many
minima, and one with a kink. Each carries its analytic gradient, a
customary starting point, its known minimizer and two flags describing
the difficulty it poses. They are exported for testing any optimizer,
not only this package's;
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
runs the whole battery.

## Usage

``` r
test_problems(which = NULL)
```

## Arguments

- which:

  A character vector of problem names, or `NULL`. `NULL`, the default,
  returns all eight in the order of the table above; a character vector
  returns just those, in the order given, so
  `test_problems(c("beale", "sphere"))` comes back beale first. A name
  that is not one of the eight raises an error listing the eight.

## Value

A named list, one element per problem, each itself a list of eight
components:

- `name`:

  character, the same as the element's name.

- `fn`:

  the objective, a function of a numeric vector of length `p` returning
  one number.

- `gr`:

  its analytic gradient, returning a vector of length `p`.

- `par`:

  numeric of length `p`, the customary starting point.

- `solution`:

  numeric of length `p`, the minimizer. For `himmelblau` this is the one
  of the four nearest the start.

- `value`:

  numeric, the minimum. Zero for all eight.

- `multimodal`:

  logical, `TRUE` for `himmelblau` and `rastrigin`.

- `smooth`:

  logical, `FALSE` for `abs_sum` alone.

## The eight

|  |  |  |  |
|----|----|----|----|
| **name** | **p** | **start** | **what it tests** |
| `sphere` | 3 | (1.3, -0.7, 0.8) | nothing; a run that fails here is broken |
| `rosenbrock` | 2 | (-1.2, 1) | a curved valley, Hessian condition 2510 at the solution |
| `booth` | 2 | (0, 0) | a well-conditioned quadratic in disguise |
| `beale` | 2 | (1, 1) | a narrow valley, Hessian condition 163 |
| `powell` | 4 | (3, -1, 0, 1) | a **singular** Hessian at the solution |
| `himmelblau` | 2 | (0, 0) | four minima, all of value zero |
| `rastrigin` | 2 | (0.4, -0.4) | 121 local minima on the usual box |
| `abs_sum` | 3 | (0, 0, 0) | a kink at the solution itself |

The starting points are the customary ones, which for `rosenbrock` and
`powell` are the deliberately unhelpful points those functions were
designed around. The two hardest are hard for opposite reasons.
`powell`'s Hessian at the origin has eigenvalues 202, 20, 0 and 0, so a
second-order method has no curvature to read in two of its four
directions and converges linearly where it usually converges
quadratically. `rastrigin` has curvature everywhere and 121 places to
stop, so what decides the answer is where the run began.

## The two flags

`multimodal` marks a surface with more than one local minimum. A local
method that reaches a different one there is behaving correctly, so
scoring it against `solution` would report a failure that is not one.
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
therefore scores every problem on the value reached and labels these two
in its `note` column.

`smooth` is `FALSE` for `abs_sum` alone. Its gradient at the solution is
exactly zero, because `sign(0)` is zero, and at every neighboring point
the subgradient has max-norm one. A method that tests \\\nabla f\\ will
therefore arrive at the answer and be unable to certify it; the
derivative-free and non-smooth methods report a measure of their own
instead, tested by
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md).

## Every solution has value zero, and that hides an effect

A line search accepts a step only when the objective decreases by a
definite amount, and near a minimum that decrease is about \\\lVert g
\rVert^{2} / (2\lambda)\\ for a curvature \\\lambda\\. Once it falls
below the rounding of the objective itself, about \\\varepsilon \lvert f
\rvert\\, no step in any direction can be verified and the search
refuses all of them. The smallest gradient a run can reach is therefore
of order

\$\$\lVert g \rVert\_{\text{floor}} \approx \sqrt{2 \lambda \varepsilon
\lvert f(x^{\ast}) \rvert},\$\$

which grows with the value at the solution and is exactly zero for every
problem here. A log-likelihood is the opposite case, of order one at its
optimum, so an optimizer that reaches `1e-15` on this battery may stop
at `1e-8` there. The defaults of
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
allow for that.

Adding a constant to any of these objectives moves neither the minimizer
nor the gradient, and reproduces the effect on demand. Conjugate
gradients on `rosenbrock`, asked for `crit_grad(1e-14)` with a budget of
20000 iterations, reaches a max-norm gradient of `4.4e-15` at \\f^{\ast}
= 0\\ and stops at `4.4e-08`, `2.8e-06` and `6.5e-05` once the constants
\\1\\, \\10^{3}\\ and \\10^{6}\\ are added. Only the first run reports
`converged = TRUE`; in the other three the tolerance asked for was never
reachable.

## See also

[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md),
which runs the battery and reports what each method reached, and
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for the rule the non-smooth problem needs.

## Examples

``` r
names(test_problems())
#> [1] "sphere"     "rosenbrock" "booth"      "beale"      "powell"    
#> [6] "himmelblau" "rastrigin"  "abs_sum"   

# Every problem's stated solution really is a stationary point of value zero.
P <- test_problems()
stopifnot(all(vapply(P, function(q) q$fn(q$solution), 0) == 0))

# BFGS solves the curved valley from the unhelpful start.
p <- test_problems("rosenbrock")[[1]]
fit <- minimize(bfgs(), p$fn, p$par, gr = p$gr)
all.equal(fit@par, p$solution, tolerance = 1e-6)
#> [1] TRUE

# Powell's difficulty is a singular Hessian at the solution: the quartic
# terms contribute nothing to the curvature at the origin, so two of the
# four eigenvalues are exactly zero and a second-order method has no
# curvature to read in those directions.
H <- matrix(c(2, 20, 0, 0,
              20, 200, 0, 0,
              0, 0, 10, -10,
              0, 0, -10, 10), 4, 4)
eigen(H, only.values = TRUE)$values
#> [1] 202  20   0   0

# abs_sum has a zero gradient at the solution and a subgradient of max-norm
# one at every neighboring point, so a gradient rule cannot certify it.
k <- test_problems("abs_sum")[[1]]
k$gr(k$solution)
#> [1] 0 0 0
k$gr(k$solution + 1e-9)
#> [1] 1 1 1

# A name that is not one of the eight is refused, and the message lists them.
try(test_problems("banana"))
#> Error : No such test problem: banana. Available: sphere, rosenbrock, booth, beale, powell, himmelblau, rastrigin, abs_sum
```
