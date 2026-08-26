# S7 Class for a Starting-Value Generator

The abstract parent of
[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
and
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md).
A starter stands in for the vector of starting values: it says how the
values are to be produced and carries how many of them there are, and
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
turns it into an actual vector before dispatch. Every optimizer
therefore accepts one, including an optimizer written outside the
package, and no method needs to know that starters exist.

## Usage

``` r
starter(npar = NULL)
```

## Arguments

- npar:

  The number of parameters, an integer, or `NULL` to have
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
  work it out from the bounds or from the objective.

## Value

An S7 object. The class is abstract, so every value is an object of one
of its subclasses.

## Details

The class is abstract and carries one property, `npar`, an integer or
`NULL`. A subclass needs a method for
[`starting_values()`](https://statmodels7.github.io/optimizers7/reference/starting_values.md)
and nothing else.
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
decides whether `par` is a vector to use as given or an object to
resolve by testing inheritance from this class, so a starter of your own
has to be parented here; a class carrying a
[`starting_values()`](https://statmodels7.github.io/optimizers7/reference/starting_values.md)
method but parented elsewhere is refused with
`'par' must be a numeric vector of starting values`.

This is the third of the package's extension points and it is exported
for the same reason as the other two,
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
for a stopping rule and
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
for an algorithm: the generic is of no use without the class it
dispatches on.

## See also

[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md),
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md),
[`starting_values()`](https://statmodels7.github.io/optimizers7/reference/starting_values.md).

## Examples

``` r
# The class is abstract, so it cannot be instantiated directly...
try(starter(npar = 3))
#> Error in new_object(S7_object(), npar = npar) : 
#>   Can't construct an object from abstract class <starter>

# ...but a subclass with a starting_values() method is a starter, and
# minimize() resolves it exactly as it resolves the two shipped ones.
Grid <- S7::new_class("Grid", parent = starter)
S7::method(starting_values, Grid) <- function(starter, npar)
  seq(-1, 1, length.out = npar)

starting_values(Grid(), 3)
#> [1] -1  0  1

f <- function(p) sum((p - 1:3)^2)
minimize(bfgs(), f, Grid(npar = 3))@par
#> [1] 1 2 3
```
