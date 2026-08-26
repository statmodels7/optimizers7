# The optimizer Class Object

Returns the
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
class itself, for the same reason
[`criterion_class()`](https://statmodels7.github.io/optimizers7/reference/criterion_class.md)
exists: the three functions that accept an arbitrary optimizer and must
check what they were given all take an argument named `optimizer`, which
shadows the class inside their body.

## Usage

``` r
optimizer_class()
```

## Value

The
[`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
class object, an `S7_class`.

## Details

The callers are
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md),
which wraps an optimizer,
[`chain()`](https://statmodels7.github.io/optimizers7/reference/chain.md),
which holds a list of them, and
[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md),
which tests one. Each rejects a non-optimizer with a message of its own;
without the check the failure would surface later, inside
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
dispatch.
