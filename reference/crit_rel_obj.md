# Stop When the Objective Stops Moving (Relative)

The rule \\\lvert f\_{new} - f\_{old} \rvert \< \texttt{tol}\\(\lvert
f\_{old} \rvert + \texttt{tol})\\.

## Usage

``` r
crit_rel_obj(tol = 1e-12)
```

## Arguments

- tol:

  Numeric tolerance, a single positive number. Defaults to `1e-12`,
  tighter than
  [`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)'s
  because the quantity compared is a ratio.

## Value

An S7 object of class
[CritRelObj](https://statmodels7.github.io/optimizers7/reference/CritRelObj-class.md),
inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

## Details

The `+ tol` in the denominator is a floor, and it is load-bearing: an
objective whose optimum sits at zero would otherwise be compared against
a vanishing scale, and the rule would either never fire or fire at once.

It was in the gradient methods' default rule until version 0.6.0 and is
not any more, because it never fired there: measured over the package's
own
[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md),
six methods on eight problems, the default with it and the default
without it agree on every flag, every evaluation count and every
reported point. It remains useful where an objective's scale is not
known in advance, which is exactly where
[`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)
is hard to set.

## See also

[`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)
for the version in the objective's own units,
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
for the default rule it left.

## Examples

``` r
crit_rel_obj()
#> <criterion> |df| < 1e-12 (relative)

# The same absolute change, at two objective scales.
c(small = crit_met(crit_rel_obj(1e-6), list(f_new = 1, f_old = 1 + 1e-7)),
  large = crit_met(crit_rel_obj(1e-6),
                   list(f_new = 1e6, f_old = 1e6 + 1e-7)))
#> small large 
#>  TRUE  TRUE 

# And the floor, which keeps an optimum at zero usable.
crit_met(crit_rel_obj(), list(f_new = 1e-30, f_old = 0))
#> [1] TRUE
```
