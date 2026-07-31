# Stop When the Objective Stops Moving (Absolute)

The rule \\\lvert f\_{new} - f\_{old} \rvert \< \texttt{tol}\\.

## Usage

``` r
crit_abs_obj(tol = 1e-10)
```

## Arguments

- tol:

  Numeric tolerance. Defaults to `1e-10`.

## Value

A
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
object.

## See also

[`crit_rel_obj`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md)

## Examples

``` r
crit_abs_obj()
#> <criterion> |df| < 1e-10
```
