# S7 Class for Line Searches

The abstract parent of the three line searches. A line search answers
one question, how far to travel along a direction another piece of the
method has already chosen, and it answers it as an object, so Newton,
BFGS, L-BFGS, conjugate gradients, gradient descent and Barzilai-Borwein
share one carefully written answer and a caller can replace it without
touching the method.

## Usage

``` r
line_search(label = character(0))
```

## Arguments

- label:

  A short character label, shown when the optimizer carrying it is
  printed.

## Value

An S7 object of class `line_search`, carrying `label`. The class is
abstract, so every value is an object of one of its three subclasses.

## Details

What every subclass guarantees is *sufficient decrease*, Armijo's
condition

\$\$f(x + \alpha d) \le f(x) + c_1 \alpha\\ g^\top d, \qquad 0 \< c_1 \<
1,\$\$

which asks for a fraction \\c_1\\ of the decrease the linear model
predicts, and so rules out steps that shrink the objective by an amount
vanishing faster than the step itself.
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
adds the strong curvature condition

\$\$\lvert \nabla f(x + \alpha d)^\top d \rvert \le c_2 \lvert g^\top d
\rvert, \qquad c_1 \< c_2 \< 1,\$\$

which excludes steps too short to have moved the directional derivative,
and is the guarantee
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
needs for its secant pair to carry usable curvature.
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
keeps Armijo's condition and replaces \\f(x)\\ by the worst of the last
few values.

The class is abstract; use
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md),
[`wolfe()`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
or
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md).
A fourth would need a branch in the compiled loop, so the set is not
extensible from outside the package as the criteria are.

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

# A method takes whichever it is given. On a problem where nothing goes
# wrong the two cost about the same; where they differ is on a method that
# needs the curvature condition, and on how a non-smooth problem is
# handled.
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

# The non-smooth problem of the battery, where Armijo is the better of the
# two and Wolfe stops without converging.
k <- test_problems("abs_sum")[[1]]
rbind(armijo = c(value = minimize(bfgs(line_search = armijo()), k$fn,
                                  k$par, gr = k$gr)@value),
      wolfe  = c(value = minimize(bfgs(line_search = wolfe()), k$fn,
                                  k$par, gr = k$gr)@value))
#>               value
#> armijo 6.837072e-08
#> wolfe  1.264625e-02
```
