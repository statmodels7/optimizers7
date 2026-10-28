# The Standard Test Problems

Returns eight standard test functions of the optimization literature: a
quadratic, a curved valley, three difficult polynomials, two functions
with several minima, and one with a kink. Each carries its analytic
gradient and Hessian, a customary starting point, its known minimizer
and two flags describing the difficulty that it poses. They are exported
for testing any optimizer, including one outside the package;
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
runs the whole battery and passes the Hessian to every method, so a
method that requires one is checked like the others.

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

A named list with one element per problem, and each element is a list of
nine components:

- `name`:

  character, the same as the element's name.

- `fn`:

  the objective, a function of a numeric vector of length `p` returning
  one number.

- `gr`:

  its analytic gradient, returning a vector of length `p`.

- `hess`:

  its analytic Hessian, returning a `p x p` matrix. For `abs_sum` it is
  the zero matrix, the Hessian wherever it exists.

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

  logical, `FALSE` only for `abs_sum`.

## The eight

|  |  |  |  |
|----|----|----|----|
| name | p | start | what it tests |
| `sphere` | 3 | (1.3, -0.7, 0.8) | no difficulty; a method that fails here is defective |
| `rosenbrock` | 2 | (-1.2, 1) | a curved valley, Hessian condition 2508 at the solution |
| `booth` | 2 | (0, 0) | a well-conditioned quadratic written as a sum of squares |
| `beale` | 2 | (1, 1) | a narrow valley, Hessian condition 162 |
| `powell` | 4 | (3, -1, 0, 1) | a singular Hessian at the solution |
| `himmelblau` | 2 | (0, 0) | four minima, all of value zero |
| `rastrigin` | 2 | (0.4, -0.4) | 121 local minima on the usual box |
| `abs_sum` | 3 | (0, 0, 0) | a kink at the solution itself |

The starting points are the customary ones, which for `rosenbrock` and
`powell` are the difficult points for which those functions were
designed. The Hessian of `powell` at the origin has eigenvalues 202, 20,
0 and 0, so a second-order method has no curvature in two of its four
directions and converges linearly where it usually converges
quadratically. `rastrigin`, by contrast, has 121 local minima on the
usual box, so the answer depends on the point at which the run began.

## The two flags

`multimodal` marks a surface with more than one local minimum. A local
method that reaches a different minimum there is behaving correctly, and
comparing its point with `solution` would report a failure where there
is none.
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
therefore scores every problem on the value reached and labels these two
in its `note` column.

`smooth` is `FALSE` only for `abs_sum`. Its gradient at the solution is
exactly zero, because `sign(0)` is zero, and at every neighboring point
the subgradient has max-norm one. A method that tests \\\nabla f\\
therefore reaches the solution without being able to certify it, unless
an iterate lands exactly on it; the derivative-free and non-smooth
methods report a measure of their own instead, tested by
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md).

## The gradient floor when every solution has value zero

A line search accepts a step only when the objective decreases by a
definite amount, and near a minimum that decrease is about \\\lVert g
\rVert^{2} / (2\lambda)\\ for a curvature \\\lambda\\. Once it falls
below the rounding error of the objective, about \\\varepsilon \lvert f
\rvert\\, the search cannot verify a trial step and rejects every one.
The smallest gradient that a run can reach is therefore of order

\$\$\lVert g \rVert\_{\text{floor}} \approx \sqrt{2 \lambda \varepsilon
\lvert f(x^{\ast}) \rvert},\$\$

which grows with the value at the solution and is zero for every problem
here. A log-likelihood divided by the number of observations is of order
one at its optimum, so an optimizer that reaches `1e-15` on this battery
may stop at `1e-8` there. The default of
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
allows for that.

Adding a constant to any of these objectives changes neither the
minimizer nor the gradient, and it reproduces the effect:
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md) on
`rosenbrock` with `crit_grad(1e-14)` and `maxit = 20000` converges at
\\f^{\ast} = 0\\, and stops without converging, at a larger gradient,
when a constant of \\1\\, \\10^{3}\\ or \\10^{6}\\ is added (see the
examples of
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)).

## See also

[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md),
which runs the battery and reports the values reached by a method, and
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for the rule of the non-smooth methods.

## Examples

``` r
names(test_problems())
#> [1] "sphere"     "rosenbrock" "booth"      "beale"      "powell"    
#> [6] "himmelblau" "rastrigin"  "abs_sum"   

# The stated solution of every problem has value zero.
P <- test_problems()
stopifnot(all(vapply(P, function(q) q$fn(q$solution), 0) == 0))

# BFGS solves the curved valley from the customary start.
p <- test_problems("rosenbrock")[[1]]
fit <- minimize(bfgs(), p$fn, p$par, gr = p$gr)
all.equal(fit@par, p$solution, tolerance = 1e-6)
#> [1] TRUE

# Powell's difficulty is a singular Hessian at the solution: the quartic
# terms contribute nothing to the curvature at the origin, so two of the
# four eigenvalues are exactly zero (they print with rounding error of
# order 1e-14) and a second-order method has no curvature to read in those
# directions.
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

# A name that is not one of the eight is rejected, and the message lists
# them.
try(test_problems("banana"))
#> Error : There is no test problem named banana. The available problems are sphere, rosenbrock, booth, beale, powell, himmelblau, rastrigin, abs_sum
```
