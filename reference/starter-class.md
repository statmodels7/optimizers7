# S7 Class for a Starting-Value Generator

This class is the abstract parent of
[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
and
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md).
A starter takes the place of the vector of starting values: it specifies
how the values are produced and how many there are, and
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
turns it into a vector before dispatch. Every optimizer therefore
accepts one, including an optimizer written outside the package, and the
methods do not need code of their own for starters.

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
tests inheritance from this class to determine whether `par` is a vector
to use as given or an object to resolve, so a new starter has to inherit
from it; a class that has a
[`starting_values()`](https://statmodels7.github.io/optimizers7/reference/starting_values.md)
method but inherits from another class is rejected with the error
`'par' must be a numeric vector of starting values`.

A starter is an extension point of the package, like
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
for a stopping rule and
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
for an algorithm, and its class is exported for the same reason: the
generic is of no use without the class on which it dispatches.

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
# minimize() resolves it as it resolves the two shipped ones.
Grid <- S7::new_class("Grid", parent = starter)
S7::method(starting_values, Grid) <- function(starter, npar)
  seq(-1, 1, length.out = npar)

starting_values(Grid(), 3)
#> [1] -1  0  1

f <- function(p) sum((p - 1:3)^2)
minimize(bfgs(), f, Grid(npar = 3))@par
#> [1] 1 2 3
```
