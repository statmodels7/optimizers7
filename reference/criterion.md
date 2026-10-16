# S7 Class for Convergence Criteria

A criterion is the abstract parent of every stopping rule. An optimizer
carries one, and a caller may replace it, combine several, or write a
new kind: a class that inherits from this one and has a method for
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
is a stopping rule, and every algorithm in the package consults it.

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

A rule is an object that implements one generic, so a rule written by
the caller is treated in the same way as a shipped one. A new rule needs
no compiled code, whereas a new algorithm or a new kind of line search
needs a branch in the compiled loop.

## The two generics

[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
returns whether the run should stop, given the state of the iteration
just completed.
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)
names the `state` components that the rule reads, so that an optimizer
unable to supply one of them rejects the rule when the run starts. The
base method for
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)
returns [`character()`](https://rdrr.io/r/base/character.html), so a
rule that reads only the objective needs no method.

## Combining

[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
and
[`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md)
combine rules and are themselves criteria, so combinations nest and
their labels nest with them.
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md),
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md)
and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
default to a disjunction of three rules.

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

# ...but a class that inherits from it is a criterion, including a rule
# written by the caller.
Tiny <- S7::new_class("Tiny", parent = criterion,
                      properties = list(tol = S7::class_numeric))
S7::method(crit_met, Tiny) <- function(criterion, state)
  state$f_new < criterion@tol
crit_met(Tiny(label = "f < 1e-6", tol = 1e-6), list(f_new = 1e-9))
#> [1] TRUE
```
