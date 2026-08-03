# Nonmonotone Backtracking

Armijo backtracking that compares against the worst of the last few
objective values rather than against the current one, so a step is
allowed to make things worse now in order to be better placed later.

## Usage

``` r
nonmonotone(c1 = 1e-04, shrink = 0.5, memory = 10, max_step = 30)
```

## Arguments

- c1:

  Sufficient-decrease constant. Defaults to `1e-4`.

- shrink:

  Factor applied to the step on each backtrack. Defaults to `0.5`.

- memory:

  How many earlier values to look back over. Defaults to `10`; `0` makes
  this ordinary
  [`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md).

- max_step:

  Most backtracks before the search gives up. Defaults to `30`.

## Value

A
[`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
object.

## Details

The condition is Grippo, Lampariello and Lucidi's: \$\$f(x_k + s d_k)
\le \max\_{0 \le j \le m} f(x\_{k-j}) + c_1 s\\ g_k^ op d_k,\$\$ which
is Armijo's with the reference replaced by the largest of the last
\\m+1\\ values. Every step it accepts improves on the worst of recent
memory; none is required to improve on the present.

### Purpose

Some methods are efficient *because* of steps that make the objective
worse. [`bb`](https://statmodels7.github.io/optimizers7/reference/bb.md)
is the clear case: its step length is a curvature estimate taken from
the last secant pair, and following that estimate faithfully means
occasionally going somewhere higher in order to be aligned with the
curvature when it matters. An Armijo condition forbids exactly those
steps and backtracks until it finds a shorter one, which is safe and is
also most of what the method was for.

The cost is that the guarantee weakens. A monotone method cannot cycle,
because the objective is a decreasing sequence bounded below; a
nonmonotone one needs the finite memory to play that role, and the
convergence result is correspondingly more delicate. Use it where a
method asks for it, not as a faster default.

### Compatibility

The curvature condition is a statement about the gradient at the trial
point and has nothing to say about which value the decrease is measured
against, so a nonmonotone Wolfe search would be a different object
rather than an option on this one. There is not one here.

## References

Grippo, L., Lampariello, F. and Lucidi, S. (1986). A nonmonotone line
search technique for Newton's method. *SIAM Journal on Numerical
Analysis* **23**, 707–716.

Raydan, M. (1997). The Barzilai and Borwein gradient method for the
large scale unconstrained minimization problem. *SIAM Journal on
Optimization* **7**, 26–33.

## See also

[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md),
[`bb`](https://statmodels7.github.io/optimizers7/reference/bb.md)

## Examples

``` r
nonmonotone()
#> <line_search> nonmonotone backtracking (memory = 10)
nonmonotone(memory = 5)
#> <line_search> nonmonotone backtracking (memory = 5)

# what it buys the method it was added for
f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))
minimize(bb(line_search = armijo()), f, c(-1.2, 1), gr = gr)@iterations
#> [1] 82
minimize(bb(), f, c(-1.2, 1), gr = gr)@iterations
#> [1] 68
```
