# Stop When the Parameters Stop Moving (Relative)

The rule \\\max_j \lvert x_j^{new} - x_j^{old}\rvert / (\lvert
x_j^{old}\rvert + \texttt{tol}) \< \texttt{tol}\\.

## Usage

``` r
crit_rel_par(tol = 1e-08)
```

## Arguments

- tol:

  Numeric tolerance, a single positive number. Defaults to `1e-8`.

## Value

An S7 object of class
[CritRelPar](https://statmodels7.github.io/optimizers7/reference/CritRelPar-class.md),
inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

## Details

The `+ tol` in each denominator is a floor, as in
[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md):
a coordinate sitting at zero would otherwise be divided by nothing. The
same number therefore serves as the floor and as the threshold.

## See also

[`crit_abs_par()`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md)
for the version in the parameters' own units.

## Examples

``` r
crit_rel_par()
#> <criterion> |dx| < 1e-08 (relative)

# The same relative change at two very different scales.
c(large = crit_met(crit_rel_par(1e-6), list(x_new = 1000,
                                            x_old = 1000.0001)),
  small = crit_met(crit_rel_par(1e-6), list(x_new = 1e-3,
                                            x_old = 1e-3 + 1e-10)))
#> large small 
#>  TRUE  TRUE 
```
