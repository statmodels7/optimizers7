# Stop When the Objective Stops Moving (Relative)

The rule \\\lvert f\_{new} - f\_{old} \rvert \< \texttt{tol}\\(\lvert
f\_{old} \rvert + \texttt{tol})\\.

## Usage

``` r
crit_rel_obj(tol = 1e-12)
```

## Arguments

- tol:

  Numeric tolerance. Defaults to `1e-12`.

## Value

A
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
object.

## Details

The `+ tol` in the denominator is a floor, and it is load-bearing: an
objective whose optimum is at zero would otherwise be compared against a
vanishing scale.

## See also

[`crit_abs_obj`](https://statmodels7.github.io/optimizers7/reference/crit_abs_obj.md)

## Examples

``` r
crit_rel_obj()
#> <criterion> |df| < 1e-12 (relative)
```
