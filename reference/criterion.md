# S7 Class for Convergence Criteria

The abstract parent of every stopping rule. An optimizer carries one,
and a caller may replace it, combine several, or write a new kind: a
class inheriting from this one with a method for
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
is a stopping rule and every algorithm in the package will consult it.

## Usage

``` r
criterion(label = character(0))
```

## Arguments

- label:

  A short character label, reported as `criterion_met` when the rule
  fires and shown when the optimizer carrying it is printed.

## Value

An S7 object of class `criterion`, carrying `label`. The class is
abstract, so every value is an object of one of its subclasses or of a
subclass written by the caller.

## Details

The alternative would be an argument taking a string and a `switch`
inside every algorithm, which fixes the set of rules at the moment the
package is written and lets nothing outside add to it. Here a rule is an
object implementing one generic, so a rule of your own is treated
exactly as a shipped one. This is the package's most open extension
point: the algorithms need a branch in compiled code and the line
searches do too, but a criterion needs neither.

## The two generics

[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
answers whether the run should stop, given the state of the iteration
just completed.
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)
names the `state` components the rule reads, so that an optimizer unable
to fill one in can refuse the rule when the run starts instead of
accepting one that never fires. The base method for
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)
returns [`character()`](https://rdrr.io/r/base/character.html), so a
rule reading only the objective needs no method at all.

## Combining

[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
and
[`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)
combine rules and are themselves criteria, so combinations nest and
their labels nest with them. The gradient methods default to a
disjunction of three.

## See also

[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md),
[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md),
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md),
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)

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
