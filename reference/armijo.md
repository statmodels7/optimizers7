# Backtracking Line Search with the Armijo Condition

Halves the step until the objective decreases by enough: \$\$f(x +
\alpha d) \leq f(x) + c_1 \alpha\\ g^\top d .\$\$ Evaluates the
objective at each trial point and never the gradient, which makes it the
cheap line search and the one a method with no curvature approximation
to protect should use.

## Usage

``` r
armijo(c1 = 1e-04, shrink = 0.5, max_step = 30, resolution = 0)
```

## Arguments

- c1:

  Sufficient-decrease constant, strictly inside \\(0, 1)\\. Defaults to
  `1e-4`, the conventional value: it demands a decrease, but only a tiny
  fraction of what the linear model predicts, so it almost never rejects
  a sensible step.

- shrink:

  Factor applied on each backtrack, strictly inside \\(0, 1)\\. Defaults
  to `0.5`.

- max_step:

  Maximum **backtracks** before the search gives up, a positive whole
  number. Defaults to 30. The name says step and the quantity is a
  count.

- resolution:

  The smallest difference in the objective that means anything, in the
  objective's own units, or a function of no arguments returning it
  where it moves as the run goes. Defaults to `0`, which does not ask
  the question. See the section below.

## Value

An S7 object of class
[ArmijoSearch](https://statmodels7.github.io/optimizers7/reference/ArmijoSearch-class.md),
inheriting from
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md),
to be passed as an optimizer's `line_search`.

## Details

The \\c_1 \alpha\\ g^\top d\\ term is essential, and dropping it to test
merely \\f\_{new} \le f\\ is a real defect. On a quadratic with unit
step the gradient update reflects the iterate through the minimum,
leaving the objective *exactly* unchanged; the weak test accepts that
step, the iterate oscillates for ever, and a stopping rule watching the
objective sees no change and reports convergence at a point that is not
a minimum.

The search evaluates the objective at trial points and never the
gradient. That is enough for a method that has only to make progress. A
quasi-Newton method needs more, and
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
supplies it.

## What the objective can resolve

An objective computed by a procedure instead of by a formula returns
slightly different values for the same argument: a fit warm-started from
wherever the last evaluation ended, a quadrature whose panels move, a
simulation. Below that spread its values carry no information, and a
search asked to verify a smaller decrease backtracks to exhaustion.

The test is made **once**, before any trial is paid for, on the
improvement the method's own linear model predicts over the full step,
\\\alpha_0 \lvert g^\top d\rvert\\. Where that is below `resolution` the
search returns immediately, and the run reports that the point is
optimal to the accuracy the objective has. That is a weaker statement
than a stopping rule being met, and it is reported in different words.

**It is not asked inside the backtracking loop, and that restriction is
what keeps it safe.** There the two situations cannot be told apart,
since \\x + \alpha d \to x\\ as the step shrinks and the objective stops
resolving the change whether the point is optimal or the direction is
wrong. Tested at the full step they separate: a bad direction predicts a
large improvement, does not obtain it, and is still reported as the
failure it is. Measured, with the test inside the loop a mis-stated
gradient at a point nowhere near stationary was promoted to a converged
run.

The quantity is the predicted decrease and not the Armijo demand \\c_1
\alpha_0 \lvert g^\top d\rvert\\, which is four orders smaller and would
fire where the method still had real progress to make.

### A resolution that moves

Where the objective settles as the run goes, a fit warm-started from the
previous evaluation locating its own answer better each time, the
resolution at the start is the reading from the worst point of the whole
run. Passing a **function** of no arguments instead of a number has it
asked again at every iteration, once per invocation of the search rather
than once per trial, so it costs one call an iteration. What the
function returns is the resolution in force for the step about to be
taken; a value that is not finite and positive is read as `0`, which
asks nothing.

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
for the shared contract.

## Examples

``` r
armijo()
#> <line_search> Armijo backtracking (c1 = 1e-04)
armijo(shrink = 0.2, max_step = 10)
#> <line_search> Armijo backtracking (c1 = 1e-04)

minimize(gd(line_search = armijo(shrink = 0.2)),
         function(p) sum((p - 1:2)^2), c(0, 0))@par
#> [1] 0.9999972 1.9999943

# The cheapness shows where the method has to backtrack often. Under
# Barzilai-Borwein, which relies on steps that go uphill, an Armijo
# condition rejects exactly those and the run costs more than twice as
# much.
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))
c(armijo      = minimize(bb(line_search = armijo()), f, c(-1.2, 1),
                         gr = g)@counts[["f"]],
  nonmonotone = minimize(bb(), f, c(-1.2, 1), gr = g)@counts[["f"]])
#>      armijo nonmonotone 
#>         154          67 

# Constants outside their intervals are refused by name.
try(armijo(c1 = 1))
#> Error : 'c1' must be a single number strictly between 0 and 1.
try(armijo(shrink = 1.5))
#> Error : 'shrink' must be a single number strictly between 0 and 1.
try(armijo(max_step = 2.5))
#> Error : 'max_step' must be a single positive whole number.
```
