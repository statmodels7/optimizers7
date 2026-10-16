# Whether an Optimizer Takes Box Bounds

Returns `TRUE` when the optimizer applies the `lower` and `upper`
arguments of
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
so that
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
tests them, and `FALSE` for a method that takes its constraint in
another way.

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

Every method in the package except
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
removes bounds by reparametrization and returns `TRUE`, which is the
default.
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
receives a constraint inside the proximal operator, where it combines
with the non-smooth term, so bounds beside the objective would be a
second, conflicting way of imposing the same constraint.

A wrapper returns the answer of what it contains.
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md)
is bounded only when every stage is, since the bounds are passed to all
of them, and
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
returns the answer of its inner optimizer.

[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
reads the declaration before its seventh check, so a method that returns
`FALSE` is not failed for rejecting bounds.

## See also

[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
for the other declaration generic,
[`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
for what a method that returns `TRUE` receives,
[`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
for the one that returns `FALSE`.

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
