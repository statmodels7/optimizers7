# S7 Class for Line Searches

The abstract parent of the three line searches. A line search decides
how far to travel along a direction that the method has already chosen.
It is an object, so the descent methods
([`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md),
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md))
share one implementation, and a caller can replace it without changing
the method.

## Usage

``` r
line_search(label = character(0))
```

## Arguments

- label:

  A short character label, shown when the optimizer that contains the
  search is printed.

## Value

An S7 object of class `line_search`, carrying `label`. The class is
abstract, so every value is an object of one of its three subclasses.

## Details

Every subclass guarantees *sufficient decrease*. For
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
and
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
this is Armijo's condition

\$\$f(x + \alpha d) \le f(x) + c_1 \alpha\\ g^\top d, \qquad 0 \< c_1 \<
1,\$\$

which requires at least the fraction \\c_1\\ of the decrease that the
linear model predicts, and so rules out steps that reduce the objective
by an amount that vanishes faster than the step itself.
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
adds the strong curvature condition

\$\$\lvert \nabla f(x + \alpha d)^\top d \rvert \le c_2 \lvert g^\top d
\rvert, \qquad c_1 \< c_2 \< 1,\$\$

which excludes steps too short to have changed the directional
derivative, and which gives
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
a secant pair with usable curvature.
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
keeps the form of Armijo's condition and replaces \\f(x)\\ by the
largest of the last few values.

The class is abstract; use
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md),
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
or
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md).
A fourth search would need a branch in the compiled loop, so, unlike the
criteria, the set of line searches cannot be extended from outside the
package.

## Notation

\\x\\ is the current point, \\g = \nabla f(x)\\ the gradient there,
\\d\\ a direction satisfying \\g^\top d \< 0\\, and \\\alpha \> 0\\ the
**step length** the search returns. The vector \\\alpha d\\ is the step
taken; \\s\\ elsewhere in this package is the secant vector \\x^{+} -
x\\ of a quasi-Newton update, a different quantity.

## See also

[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md),
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md),
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)

## Examples

``` r
# Abstract: use one of the constructors.
try(line_search(label = "mine"))
#> Error in new_object(S7_object(), label = label) : 
#>   Can't construct an object from abstract class <line_search>

armijo()
#> <line_search> Armijo backtracking (c1 = 1e-04)
wolfe()
#> <line_search> strong Wolfe (c1 = 1e-04, c2 = 0.9)

nonmonotone()
#> <line_search> nonmonotone backtracking (memory = 10)

# A method uses the search that it is given. On a smooth problem the two
# cost about the same; they differ on a method that needs the curvature
# condition and on a non-smooth problem.
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))
rbind(armijo = unlist(minimize(bfgs(line_search = armijo()), f,
                               c(-1.2, 1), gr = g)@counts),
      wolfe  = unlist(minimize(bfgs(line_search = wolfe()), f,
                               c(-1.2, 1), gr = g)@counts))
#>         f  g h
#> armijo 49 40 0
#> wolfe  49 39 0

# The non-smooth problem of the battery, on which BFGS converges with
# Armijo and stops without converging with Wolfe.
k <- test_problems("abs_sum")[[1]]
rbind(armijo = c(value = minimize(bfgs(line_search = armijo()), k$fn,
                                  k$par, gr = k$gr)@value),
      wolfe  = c(value = minimize(bfgs(line_search = wolfe()), k$fn,
                                  k$par, gr = k$gr)@value))
#>               value
#> armijo 6.837072e-08
#> wolfe  1.264625e-02
```
