# Whether an Optimizer Takes Box Bounds

`TRUE` when the optimizer honours the `lower` and `upper` arguments of
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
so that
[`check_optimizer`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
tests them, and `FALSE` for a method that takes its constraint another
way.

## Usage

``` r
optimizer_bounded(optimizer)
```

## Arguments

- optimizer:

  An
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Value

A single logical.

## Details

Every method of the package removes bounds by reparametrization and
answers `TRUE`, which is the default. A proximal method is the
exception: a constraint reaches it inside the proximal operator, where
it composes with the term already there, so bounds beside the objective
would be a second and conflicting route to the same thing.

## See also

[`optimizer_provides`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)

## Examples

``` r
optimizer_bounded(bfgs())
#> [1] TRUE
```
