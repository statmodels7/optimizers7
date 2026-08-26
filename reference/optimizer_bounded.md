# Whether an Optimizer Takes Box Bounds

`TRUE` when the optimizer honors the `lower` and `upper` arguments of
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
so that
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
tests them, and `FALSE` for a method that takes its constraint another
way.

## Usage

``` r
optimizer_bounded(optimizer)
```

## Arguments

- optimizer:

  An
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

## Value

A single logical.

## Details

Every method here removes bounds by reparametrization and answers
`TRUE`, which is the default.
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
is the one exception: a constraint reaches it inside the proximal
operator, where it composes with the term already there, so bounds
beside the objective would be a second and conflicting route to the same
thing.

A wrapper answers for what it contains.
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
is bounded only when **every** stage is, the bounds being passed to all
of them.

The declaration is read by
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
before its seventh check, so a method answering `FALSE` is not failed
for refusing a box it never promised.

## See also

[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
for the other declaration generic,
[`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
for what a method that answers `TRUE` receives,
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
for the one that answers `FALSE`.

## Examples

``` r
optimizer_bounded(bfgs())
#> [1] TRUE

pg <- prox_grad(prox = function(v, t) v, g = function(b) 0)
optimizer_bounded(pg)
#> [1] FALSE

# A chain is bounded only if every stage is.
c(optimizer_bounded(chain(bfgs(), newton())),
  optimizer_bounded(chain(bfgs(), pg)))
#> [1]  TRUE FALSE
```
