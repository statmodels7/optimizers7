# Stop When the Parameters Stop Moving (Absolute)

Builds the rule \\\max_j \lvert x_j^{new} - x_j^{old} \rvert \<
\texttt{tol}\\, so the tolerance is in the units of the parameters. It
reads only the iterates, and it is one of the three terms of the default
rule of the gradient methods.

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

With box constraints the comparison is on the unconstrained scale, on
which the optimizer moves. Near a bound a short distance in \\\theta\\
corresponds to a long distance in \\\eta\\, so the rule fires later
there than the change on the parameter scale suggests.

A method whose iterate does not move at an iteration, such as
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
after a trial point that does not improve on the best one, reports a
change of zero, so the rule fires at that iteration.

## See also

[`crit_rel_par()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_par.md)
for the version relative to each coordinate,
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
for the default rule that contains this one.

## Examples

``` r
crit_abs_par()
#> <criterion> |dx| < 1e-08
crit_met(crit_abs_par(1e-6), list(x_new = c(1, 2), x_old = c(1, 2 + 1e-9)))
#> [1] TRUE
```
