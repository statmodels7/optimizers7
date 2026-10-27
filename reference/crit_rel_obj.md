# Stop When the Objective Stops Moving (Relative)

The rule is \\\lvert f\_{new} - f\_{old} \rvert \< \texttt{tol}\\(\lvert
f\_{old} \rvert + \texttt{tol})\\.

## Usage

``` r
crit_rel_obj(tol = 1e-12)
```

## Arguments

- tol:

  Numeric tolerance, a single positive number. Defaults to `1e-12`,
  smaller than the default of
  [`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)
  because the quantity compared is a ratio.

## Value

An S7 object of class
[CritRelObj](https://statmodels7.github.io/optimizers7/reference/CritRelObj-class.md),
inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

## Details

The term `+ tol` is a floor on the scale. Without it the threshold would
be zero at an objective value of zero, and since the comparison is
strict the rule could never fire at an optimum of zero.

The rule is useful where the scale of the objective is not known in
advance, which is where
[`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)
is hard to set.

## See also

[`crit_abs_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)
for the version in the units of the objective,
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
for the default rule of the gradient methods.

## Examples

``` r
crit_rel_obj()
#> <criterion> |df| < 1e-12 (relative)

# The same absolute change is tested at two objective scales.
c(small = crit_met(crit_rel_obj(1e-6), list(f_new = 1, f_old = 1 + 1e-7)),
  large = crit_met(crit_rel_obj(1e-6),
                   list(f_new = 1e6, f_old = 1e6 + 1e-7)))
#> small large 
#>  TRUE  TRUE 

# With the floor, the rule can fire at an optimum of zero.
crit_met(crit_rel_obj(), list(f_new = 1e-30, f_old = 0))
#> [1] TRUE
```
