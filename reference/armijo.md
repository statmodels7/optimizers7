# Backtracking Line Search with the Armijo Condition

Multiplies the step by `shrink` until the objective decreases by enough:
\$\$f(x + \alpha d) \leq f(x) + c_1 \alpha\\ g^\top d .\$\$ It evaluates
the objective at each trial point and never the gradient, which makes it
the cheaper line search, suited to a method that keeps no curvature
approximation.

## Usage

``` r
armijo(c1 = 1e-04, shrink = 0.5, max_step = 30, resolution = 0)
```

## Arguments

- c1:

  Sufficient-decrease constant, strictly inside \\(0, 1)\\. Defaults to
  `1e-4`, the conventional value: it requires a decrease, but only a
  small fraction of the decrease that the linear model predicts.

- shrink:

  Factor applied on each backtrack, strictly inside \\(0, 1)\\. Defaults
  to `0.5`.

- max_step:

  Maximum number of trial points, including the full step, before the
  search returns without a step; a positive whole number. Defaults
  to 30. The name suggests a length, but the quantity is a count.

- resolution:

  The smallest difference in the objective that the objective can
  resolve, in the units of the objective, or a function of no arguments
  that returns it when it changes during the run. Defaults to `0`, which
  disables the test. See the section below.

## Value

An S7 object of class
[ArmijoSearch](https://statmodels7.github.io/optimizers7/reference/ArmijoSearch-class.md),
inheriting from
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md),
to be passed as an optimizer's `line_search`.

## Details

The term \\c_1 \alpha\\ g^\top d\\ is required. A test that requires
only \\f\_{new} \le f\\ accepts steps that leave the objective
unchanged: on \\f(x) = \lVert x \rVert^2\\, whose Hessian is \\2I\\, a
unit step along \\-g\\ reflects the iterate through the minimum and
leaves the objective *exactly* unchanged. The weaker test accepts that
step, the iterate oscillates indefinitely, and a stopping rule that
reads the objective sees no change and reports convergence at a point
that is not a minimum.

A quasi-Newton method also needs the curvature condition, which
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
enforces.

## What the objective can resolve

An objective computed by a procedure, and not by a formula, can return
slightly different values for the same argument: a fit warm-started from
the point where the last evaluation ended, a quadrature whose panels
move, a simulation. Differences below that spread carry no information,
and a search that has to verify a smaller decrease backtracks until it
has used `max_step` trial points.

The test is made once, before the first trial point, on the decrease
that the linear model of the method predicts over the full step,
\\\alpha_0 \lvert g^\top d\rvert\\. When that decrease is below
`resolution` the search returns at once, and the run reports that the
point is optimal to the accuracy of the objective. That is a weaker
statement than the one made when a stopping rule is met, and the result
reports it with a different label.

The test is not applied inside the backtracking loop. There the two
situations cannot be distinguished, since \\x + \alpha d \to x\\ as the
step shrinks and the objective stops resolving the change whether the
point is optimal or the direction is wrong. At the full step they
differ: a wrong direction predicts a large decrease, does not obtain it,
and the search reports a failure. Inside the loop the test would report
an incorrect gradient at a point far from stationarity as convergence.

The quantity is the predicted decrease and not the Armijo threshold
\\c_1 \alpha_0 \lvert g^\top d\rvert\\, which at the default \\c_1 =
10^{-4}\\ is four orders of magnitude smaller and would fire while the
method can still make progress.

### A resolution that changes during the run

When the objective becomes more precise as the run proceeds (for example
a fit warm-started from the previous evaluation), the resolution
measured at the start is the worst value of the whole run. A function of
no arguments passed instead of a number is called at every iteration,
once per invocation of the search and not once per trial, so it costs
one call per iteration. Its value is the resolution for the step about
to be taken; a value that is not finite and positive is read as `0`,
which disables the test.

## References

Armijo, L. (1966). Minimization of functions having Lipschitz continuous
first partial derivatives. *Pacific Journal of Mathematics* **16**, 1–3.

Zoutendijk, G. (1970). Nonlinear programming, computational methods. In
J. Abadie (ed.), *Integer and Nonlinear Programming*, 37–86.
North-Holland, Amsterdam.

## See also

[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
for the curvature condition,
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
for the relaxed reference value,
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
for the parent class.

## Examples

``` r
armijo()
#> <line_search> Armijo backtracking (c1 = 1e-04)
armijo(shrink = 0.2, max_step = 10)
#> <line_search> Armijo backtracking (c1 = 1e-04)

minimize(gd(line_search = armijo(shrink = 0.2)),
         function(p) sum((p - 1:2)^2), c(0, 0))@par
#> [1] 0.9999972 1.9999943

# Barzilai-Borwein relies on steps that go uphill. An Armijo condition
# rejects every such step, and the run costs more than twice as many
# evaluations as under the default nonmonotone search.
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))
c(armijo      = minimize(bb(line_search = armijo()), f, c(-1.2, 1),
                         gr = g)@counts[["f"]],
  nonmonotone = minimize(bb(), f, c(-1.2, 1), gr = g)@counts[["f"]])
#>      armijo nonmonotone 
#>         154          67 

# Constants outside their intervals are rejected, and the error names
# the argument.
try(armijo(c1 = 1))
#> Error : 'c1' must be a single number strictly between 0 and 1.
try(armijo(shrink = 1.5))
#> Error : 'shrink' must be a single number strictly between 0 and 1.
try(armijo(max_step = 2.5))
#> Error : 'max_step' must be a single positive whole number.
```
