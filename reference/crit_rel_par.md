# Stop When the Parameters Stop Moving (Relative)

The rule \\\max_j \lvert x_j^{new} - x_j^{old}\rvert / (\lvert
x_j^{old}\rvert + \texttt{tol}) \< \texttt{tol}\\.

## Usage

``` r
crit_rel_par(tol = 1e-08)
```

## Arguments

- tol:

  Numeric tolerance. Defaults to `1e-8`.

## Value

A
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
object.

## See also

[`crit_abs_par`](https://statmodels7.github.io/optimizers7/reference/crit_abs_par.md)

## Examples

``` r
crit_rel_par()
#> <criterion> |dx| < 1e-08 (relative)
```
