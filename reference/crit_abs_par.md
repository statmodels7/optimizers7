# Stop When the Parameters Stop Moving (Absolute)

Builds the rule \\\max_j \lvert x_j^{new} - x_j^{old} \rvert \<
\texttt{tol}\\, so the tolerance is in the parameters' own units. It
asks nothing of the method, and it is one of the three terms in the
gradient methods' default rule.

## Usage

``` r
crit_abs_par(tol = 1e-08)
```

## Arguments

- tol:

  Numeric tolerance, a single positive number. Defaults to `1e-8`.

## Value

An S7 object of class
[CritAbsPar](https://statmodels7.github.io/optimizers7/reference/CritAbsPar-class.md),
inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

## Details

With box constraints the comparison is on the **unconstrained** scale,
that being where the optimizer moves. A coordinate approaching a bound
travels a long way in \\\eta\\ for a short way in \\\theta\\, so this
rule fires later there than a reader of the parameter scale would
expect.

## See also

[`crit_rel_par()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_par.md)
for the version scaled by each coordinate,
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
for the default rule that holds this one.

## Examples

``` r
crit_abs_par()
#> <criterion> |dx| < 1e-08
crit_met(crit_abs_par(1e-6), list(x_new = c(1, 2), x_old = c(1, 2 + 1e-9)))
#> [1] TRUE
```
