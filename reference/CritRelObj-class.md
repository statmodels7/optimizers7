# S7 Class for the Relative Objective Criterion

The rule
[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
builds. It fires when \\\lvert f\_{new} - f\_{old}\rvert \<
\texttt{tol}\\ (\lvert f\_{old}\rvert + \texttt{tol})\\, so the
comparison is against the objective's own scale and one tolerance serves
whatever units the problem is in. It declares nothing through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md).

## Usage

``` r
CritRelObj(label = character(0), tol = integer(0))
```

## Arguments

- tol:

  The tolerance, a single positive number.

## Value

An S7 object of class `CritRelObj` inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
carrying `label` and `tol`.

## Details

[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
returns `FALSE` when `f_old` is `NULL` or not finite, as at the first
iteration.

The `+ tol` in the denominator is a floor and is load-bearing: an
objective whose optimum sits at zero would otherwise be compared against
a vanishing scale, and the rule would either never fire or fire at once.

## See also

[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)
for the constructor,
[`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)
for the version in the objective's own units.

## Examples

``` r
crit_met(crit_rel_obj(1e-6), list(f_new = 1.0000001, f_old = 1.0000002))
#> [1] TRUE

# The floor is what keeps an optimum at zero usable.
crit_met(crit_rel_obj(), list(f_new = 1e-30, f_old = 0))
#> [1] TRUE
```
