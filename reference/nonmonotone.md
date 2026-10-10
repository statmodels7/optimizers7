# Nonmonotone Backtracking

Armijo backtracking that compares the trial value with the largest of
the last few objective values and not with the current one, so a step
may increase the objective.

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

  The number of earlier values, besides the current one, over which the
  reference is taken. Defaults to `10`; `0` makes the search identical
  to
  [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md).

- max_step:

  Maximum number of trial points, the full step included, before the
  search gives up. Defaults to `30`.

- resolution:

  The smallest difference in the objective that the objective can
  resolve, in the units of the objective, or a function of no arguments
  that returns it when it changes during the run. Defaults to `0`, which
  disables the test; see
  [`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
  for the purpose of the test and for the reason why it is applied at
  the full step, before the search begins.

## Value

A
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
object.

## Details

The condition is that of Grippo, Lampariello and Lucidi: \$\$f(x_k +
\alpha d_k) \le \max\_{0 \le j \le m} f(x\_{k-j}) + c_1 \alpha\\
g_k^\top d_k,\$\$ which is Armijo's condition with the reference
replaced by the largest of the last \\m+1\\ values. An accepted step
improves on the largest of the recent values, but it need not improve on
the current one.

## Purpose

Some methods are efficient *because* of steps that increase the
objective.
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) is
the main case: its step length is a curvature estimate taken from the
last secant pair, and following that estimate sometimes leads to a
higher value. An Armijo condition rejects every step that increases the
objective and backtracks until it finds a shorter one, which removes the
steps on which the method relies. On the Rosenbrock function
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md)
needs fewer than half as many evaluations of the objective under
`nonmonotone()` as under
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
(see the examples).

The cost is a weaker guarantee. A monotone method cannot cycle, since
the objective is a decreasing sequence bounded below; a nonmonotone
method relies on the finite memory for that, and its convergence result
is more delicate. The search is intended for methods such as
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md),
whose steps are not required to decrease the objective.

## memory = 0

At `memory = 0` the reference is the current value and the condition is
exactly Armijo's, so `nonmonotone(memory = 0)` and
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
produce the same run, and a comparison between the two settings isolates
the effect of the memory.

## No nonmonotone Wolfe search

The curvature condition concerns the gradient at the trial point and
does not involve the value against which the decrease is measured. A
nonmonotone Wolfe search would therefore be a fourth object and not an
option of this one, and the package does not provide one.

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
the method that uses this one by default,
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
for the parent class.

## Examples

``` r
nonmonotone()
#> <line_search> nonmonotone backtracking (memory = 10)
nonmonotone(memory = 5)
#> <line_search> nonmonotone backtracking (memory = 5)

# The effect on bb(), the method that uses it by default.
f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))
evals <- function(ls) minimize(bb(line_search = ls), f, c(-1.2, 1),
                               gr = gr)@counts[["f"]]
c(nonmonotone = evals(nonmonotone()), armijo = evals(armijo()))
#> nonmonotone      armijo 
#>          67         154 

# With memory = 0 the search is armijo(), so the comparison above
# isolates the memory.
identical(evals(nonmonotone(memory = 0)), evals(armijo()))
#> [1] TRUE
```
