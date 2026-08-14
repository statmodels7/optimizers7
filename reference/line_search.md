# S7 Class for Line Searches

How far to step along a direction, as an object.

## Usage

``` r
line_search(label = character(0))
```

## Arguments

- label:

  A short character label, used when printing.

## Value

An S7 object of class `line_search`.

## Details

A descent method produces a *direction*; how far to travel along it is a
separate question with its own theory. Separating the two is what lets
Newton, BFGS and L-BFGS share one carefully written answer instead of
each carrying its own, and what lets a user change the answer without
touching the method.

At a point \\x\\ with gradient \\g\\ and a direction \\d\\ satisfying
\\g^\top d \< 0\\, a line search returns a step length \\s \> 0\\. What
every subclass guarantees is *sufficient decrease*, Armijo's condition

\$\$f(x + s d) \le f(x) + c_1 s\\ g^\top d, \qquad 0 \< c_1 \< 1,\$\$

which asks for a fraction \\c_1\\ of the decrease the linear model
predicts and so rules out steps that shrink the objective by an amount
vanishing faster than the step itself.
[`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)
adds the strong curvature condition

\$\$\lvert \nabla f(x + s d)^\top d \rvert \le c_2 \lvert g^\top d
\rvert, \qquad c_1 \< c_2 \< 1,\$\$

which excludes steps too short to have moved the directional derivative
and is what makes the secant pair of
[`bfgs`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
carry usable curvature.

The class is abstract; use
[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
or
[`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md).

## See also

[`armijo`](https://statmodels7.github.io/optimizers7/reference/armijo.md),
[`wolfe`](https://statmodels7.github.io/optimizers7/reference/wolfe.md)

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

# A method takes whichever it is given, and the choice shows in the result.
f <- function(p) sum((p - c(1, 2))^2)
g <- function(p) 2 * (p - c(1, 2))
minimize(bfgs(line_search = armijo()), f, c(0, 0), gr = g)@counts
#> f g h 
#> 3 3 0 
minimize(bfgs(line_search = wolfe()), f, c(0, 0), gr = g)@counts
#> f g h 
#> 3 3 0 
```
