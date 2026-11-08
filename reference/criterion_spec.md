# The Compiled Form of a Stopping Rule

Returns a description of a stopping rule as nested lists, from which the
compiled descent loop builds the same rule in C++, or `NULL` where the
rule or one of its members has no compiled form.

## Usage

``` r
criterion_spec(criterion)
```

## Arguments

- criterion:

  A
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object.

## Value

A list, or `NULL`.

## Details

The rules with a compiled form are those built by
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md),
[`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md),
[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md),
[`crit_abs_par()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md),
[`crit_rel_par()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_par.md),
[`crit_never()`](https://statmodels7.github.io/optimizers7/reference/crit_never.md),
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
and
[`crit_all()`](https://statmodels7.github.io/optimizers7/reference/crit_all.md).
The class is read exactly, so that a subclass with a method of its own
for
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
is evaluated in R. Each element of the result has `type` (`"grad"`,
`"abs_obj"`, `"rel_obj"`, `"abs_par"`, `"rel_par"`, `"never"`, `"any"`
or `"all"`), and `tol`, `norm` or `criteria` where the rule has them.
The compiled rules compute the same quantities as the methods of
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md),
with the same arithmetic, and combine missing values as
[`any()`](https://rdrr.io/r/base/any.html) and
[`all()`](https://rdrr.io/r/base/all.html) do. The descent methods
([`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md),
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md),
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md))
ask the compiled form where one exists;
`options(optimizers7.compiled_rules = FALSE)` makes them ask
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
in R instead. Another package that runs the descent loop through the
headers in `inst/include/optimizers7` passes this list to
`make_compiled_rule()` in `rules.h`.

## See also

[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md),
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)

## Examples

``` r
str(criterion_spec(crit_any(crit_grad(), crit_abs_obj(1e-8))))
#> List of 2
#>  $ type    : chr "any"
#>  $ criteria:List of 2
#>   ..$ :List of 3
#>   .. ..$ type: chr "grad"
#>   .. ..$ tol : num 1e-06
#>   .. ..$ norm: chr "max"
#>   ..$ :List of 2
#>   .. ..$ type: chr "abs_obj"
#>   .. ..$ tol : num 1e-08

Tiny <- S7::new_class("Tiny", parent = criterion,
                      properties = list(tol = S7::class_numeric))
S7::method(crit_met, Tiny) <- function(criterion, state)
  abs(state$f_new) < criterion@tol
#> Overwriting method crit_met(<Tiny>)
criterion_spec(Tiny(label = "tiny", tol = 1e-3))
#> NULL
```
