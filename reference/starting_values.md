# Produce a Vector of Starting Values

Turns a starter into an actual numeric vector of length `npar`, on the
**unconstrained** scale.

## Usage

``` r
starting_values(starter, npar)
```

## Arguments

- starter:

  A
  [`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
  or
  [`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md)
  object, or a user-defined starter.

- npar:

  The number of parameters wanted.

## Value

A numeric vector of length `npar`, on the unconstrained scale.

## Details

The unconstrained scale is the one the optimizer works on when there are
bounds, and that is where a starter is entitled to be simple: zero means
the middle of an interval, one for a variance, one half for a
probability, and no value can fall outside a bound.
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
maps the result back through
[`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)
before any method sees it.

`npar` is passed by
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
which has already settled it from the starter's own `npar`, from the
length of the bounds, or by probing the objective with
[`infer_npar()`](https://statmodels7.github.io/optimizers7/reference/infer_npar.md).
A method for this generic therefore reads the argument, not
`starter@npar`.

The two shipped starters are
[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
and
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md).
A starter of a third kind is a subclass of the abstract
[starter](https://statmodels7.github.io/optimizers7/reference/starter-class.md)
class with a method for this generic and nothing else, and
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
then accepts it as `par` like either shipped one. The parent matters:
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
accepts as `par` a numeric vector or an object inheriting from
[starter](https://statmodels7.github.io/optimizers7/reference/starter-class.md),
so a class carrying a method for this generic but parented elsewhere is
refused with `'par' must be a numeric vector of starting values`.

## See also

[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
and
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md)
for the two shipped starters,
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
for where the resolution happens,
[`infer_npar()`](https://statmodels7.github.io/optimizers7/reference/infer_npar.md)
for how `npar` is worked out when nothing declares it.

## Examples

``` r
starting_values(start_zeros(), 3)
#> [1] 0 0 0

set.seed(1)
starting_values(start_runif(-2, 2), 3)
#> [1] -0.9379653 -0.5115044  0.2914135

# Resolving by hand gives exactly the vector minimize() would have used,
# so passing a starter and passing its values are the same run.
f <- function(p) sum((p - 1:3)^2)
set.seed(1); a <- minimize(bfgs(), f, start_runif(-2, 2, npar = 3))
set.seed(1); b <- minimize(bfgs(), f, starting_values(start_runif(-2, 2), 3))
identical(a@par, b@par)
#> [1] TRUE
```
