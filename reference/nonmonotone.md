# Nonmonotone Backtracking

Armijo backtracking that compares against the worst of the last few
objective values rather than against the current one, so a step is
allowed to make things worse now in order to be better placed later.

## Usage

``` r
nonmonotone(
  c1 = 1e-04,
  shrink = 0.5,
  memory = 10,
  max_step = 30,
  resolution = 0
)
```

## Arguments

- c1:

  Sufficient-decrease constant. Defaults to `1e-4`.

- shrink:

  Factor applied to the step on each backtrack. Defaults to `0.5`.

- memory:

  How many earlier values to look back over. Defaults to `10`; `0` makes
  this ordinary
  [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md).

- max_step:

  Most backtracks before the search gives up. Defaults to `30`.

- resolution:

  The smallest difference in the objective that means anything, in the
  objective's own units, or a function of no arguments returning it
  where it moves as the run goes. Defaults to `0`, which does not ask
  the question; see
  [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
  for what it is for and why the question is put at the full step,
  before the search begins.

## Value

A
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
object.

## Details

The condition is Grippo, Lampariello and Lucidi's: \$\$f(x_k + s d_k)
\le \max\_{0 \le j \le m} f(x\_{k-j}) + c_1 s\\ g_k^\top d_k,\$\$ which
is Armijo's with the reference replaced by the largest of the last
\\m+1\\ values. Every step it accepts improves on the worst of recent
memory; none is required to improve on the present.

## Purpose

Some methods are efficient *because* of steps that make the objective
worse.
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) is
the clear case: its step length is a curvature estimate taken from the
last secant pair, and following that estimate faithfully means
occasionally going somewhere higher in order to be aligned with the
curvature when it matters. An Armijo condition forbids exactly those
steps and backtracks until it finds a shorter one, which is safe and is
also most of what the method was for.

Measured on Rosenbrock,
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md)
takes 58 iterations and 67 objective evaluations under `nonmonotone()`
and 72 iterations and 154 evaluations under
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md).

The cost is that the guarantee weakens. A monotone method cannot cycle,
the objective being a decreasing sequence bounded below; a nonmonotone
one needs the finite memory to play that role, and the convergence
result is correspondingly more delicate. Use it where a method asks for
it.

## memory = 0

At `memory = 0` the reference is the current value and the condition is
Armijo's exactly. Measured, `nonmonotone(memory = 0)` and
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
give the identical run, 72 iterations and 154 evaluations on the same
problem, so a comparison between the two settings is a comparison of the
memory alone.

## No nonmonotone Wolfe

The curvature condition is a statement about the gradient at the trial
point and says nothing about which value the decrease is measured
against, so a nonmonotone Wolfe search would be a fourth object rather
than an option on this one. There is not one here.

## References

Grippo, L., Lampariello, F. and Lucidi, S. (1986). A nonmonotone line
search technique for Newton's method. *SIAM Journal on Numerical
Analysis* **23**, 707–716.

Raydan, M. (1997). The Barzilai and Borwein gradient method for the
large scale unconstrained minimization problem. *SIAM Journal on
Optimization* **7**, 26–33.

## See also

[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
for the monotone version,
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) for
the method that needs this one,
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
for the shared contract.

## Examples

``` r
nonmonotone()
#> <line_search> nonmonotone backtracking (memory = 10)
nonmonotone(memory = 5)
#> <line_search> nonmonotone backtracking (memory = 5)

# What it buys the method it was added for.
f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))
evals <- function(ls) minimize(bb(line_search = ls), f, c(-1.2, 1),
                               gr = gr)@counts[["f"]]
c(nonmonotone = evals(nonmonotone()), armijo = evals(armijo()))
#> nonmonotone      armijo 
#>          67         154 

# And that memory = 0 is armijo(), which is what makes the comparison one
# of the memory alone.
identical(evals(nonmonotone(memory = 0)), evals(armijo()))
#> [1] TRUE
```
