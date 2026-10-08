# Brent's Method in One Dimension

Minimizes a function of one variable without derivatives: a bracket of
the minimum is found by golden expansion from the starting value, and
the bracket is then reduced by Brent's combination of parabolic
interpolation and golden-section steps.

## Usage

``` r
brent(
  criterion = crit_stationary(1e-05),
  step = 1,
  max_expand = 50,
  maxit = 100,
  max_eval = Inf,
  verbose = FALSE,
  refresh = 1,
  keep_trace = FALSE
)
```

## Arguments

- criterion:

  The stopping rule, a
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object. Defaults to `crit_stationary(1e-5)`, read on the half-width of
  the bracket, so the tolerance is on the parameter's scale. A gradient
  rule is refused when the run starts.

- step:

  The length of the first trial step of the bracketing phase, a single
  positive number, default `1`.

- max_expand:

  How many golden expansions the bracketing phase may take, a single
  positive number, default `50`. A run that has not bracketed the
  minimum by then continues from the best point it has.

- maxit:

  Maximum iterations of the reduction phase. Defaults to 100.

- max_eval:

  Maximum objective evaluations. Defaults to `Inf`.

- verbose:

  Report progress? Defaults to `FALSE`.

- refresh:

  Report every this many iterations. Defaults to 1.

- keep_trace:

  Store the iteration path? Defaults to `FALSE`.

## Value

An S7 object of class
[Brent](https://statmodels7.github.io/optimizers7/reference/Brent-class.md),
inheriting from
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md),
to be handed to
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Details

The bracketing phase evaluates the starting value \\x_0\\ and \\x_0 +
h\\, turns round if the second is the higher, and steps on in the
descending direction, each step \\1.618\\ times the last, until the
value rises. The three points then satisfy \\f(b) \< f(a)\\ and \\f(b)
\< f(c)\\ with \\b\\ between \\a\\ and \\c\\, so a minimum lies inside
\\(a, c)\\.

The reduction phase is the algorithm of Brent (1973, chapter 5). At each
iteration the parabola through the three best points is tried first; its
vertex is accepted when it lies inside the bracket and the step is less
than half the step before last, and a golden-section step into the
larger part of the bracket is taken otherwise. Parabolic steps converge
superlinearly on a smooth minimum, and the golden steps bound the number
of evaluations when the function is not well modeled by a parabola. The
stationarity measure is the half-width of the bracket.

The method is for one parameter only; a starting value of any other
length is rejected. Bounds are removed by the same reparametrization the
other methods use (see
[`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)).

## References

Brent, R. P. (1973). *Algorithms for Minimization without Derivatives*.
Prentice-Hall, Englewood Cliffs.

## See also

[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
for the derivative-free method in several dimensions,
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for the rule this method reads.

## Examples

``` r
brent()
#> <optimizer> brent
#>   stop when : stationarity < 1e-05
#>   budgets   : maxit 100, evaluations Inf
#>   settings  : step = 1, max_expand = 50

# A smooth function of one variable.
f <- function(x) (x - 2)^2 + exp(-x)
r <- minimize(brent(), f, par = 0)
c(par = r@par, evaluations = r@counts[["f"]])
#>         par evaluations 
#>    2.063504   12.000000 

# Against nelder_mead() on the same function.
minimize(nelder_mead(), f, par = 0)@counts[["f"]]
#> [1] 42
```
