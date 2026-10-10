# Validate a Line Search

Checks that the value inherits from the abstract
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
class. The message names
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
as an example, since passing the string `"armijo"` instead of the object
is a common mistake.

## Usage

``` r
check_line_search(x)
```

## Arguments

- x:

  The value supplied.

## Value

Invisibly `TRUE`. Raises an error naming `line_search` otherwise.
