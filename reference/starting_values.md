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
  [`start_zeros`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
  or
  [`start_runif`](https://statmodels7.github.io/optimizers7/reference/start_runif.md)
  object, or a user-defined starter.

- npar:

  The number of parameters wanted.

## Value

A numeric vector of length `npar`.

## Details

The unconstrained scale is the one the optimizer actually works on when
there are bounds, so this is where a starter is entitled to be simple:
zero means the middle of an interval, one for a variance, one half for a
probability, and there is no way for it to fall outside a bound.
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
maps the result back through
[`bounded_transform`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)
before any method sees it.

A user-defined starter is a subclass of `starter` with a method for this
generic; nothing else is required.

## See also

[`start_zeros`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md),
[`start_runif`](https://statmodels7.github.io/optimizers7/reference/start_runif.md),
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

## Examples

``` r
starting_values(start_zeros(), 3)
#> [1] 0 0 0

set.seed(1)
starting_values(start_runif(-2, 2), 3)
#> [1] -0.9379653 -0.5115044  0.2914135
```
