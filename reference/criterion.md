# S7 Class for Convergence Criteria

A stopping rule, as an object. Every optimizer carries one, and the user
may replace it, combine several, or write a new kind.

## Usage

``` r
criterion(label = character(0))
```

## Arguments

- label:

  A short character label, used when reporting which rule fired.

## Value

An S7 object of class `criterion`. The class is abstract: use one of the
constructors, or from a user-defined subclass.

## Details

The alternative would be an argument taking a string, and a `switch`
inside every algorithm. That is exactly the arrangement this toolkit
exists to replace: it fixes the set of rules at the moment the package
is written, and nothing outside can add to it. A criterion here is an
object implementing one generic,
[`crit_met`](https://statmodels7.github.io/optimizers7/reference/crit_met.md),
so a user-defined rule is treated like a shipped one.

Criteria are combined with
[`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
and
[`crit_all`](https://statmodels7.github.io/optimizers7/reference/crit_all.md),
which are themselves criteria, so combinations nest.

## See also

[`crit_grad`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md),
[`crit_rel_obj`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md),
[`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md),
[`crit_met`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)

## Examples

``` r
# The class is abstract, so it cannot be instantiated directly...
try(criterion(label = "mine"))
#> Error in new_object(S7_object(), label = label) : 
#>   Can't construct an object from abstract class <criterion>

# ...but anything inheriting from it is a criterion, including a rule the
# package never anticipated.
Tiny <- S7::new_class("Tiny", parent = criterion,
                      properties = list(tol = S7::class_numeric))
S7::method(crit_met, Tiny) <- function(criterion, state)
  state$f_new < criterion@tol
crit_met(Tiny(label = "f < 1e-6", tol = 1e-6), list(f_new = 1e-9))
#> [1] TRUE
```
