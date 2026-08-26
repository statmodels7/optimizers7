# The criterion Class Object

Returns the
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
class itself, for code that has to test whether a value inherits from
it. The indirection exists because the functions asking that question
take an argument named `criterion`, which shadows the class inside their
own body.

## Usage

``` r
criterion_class()
```

## Value

The
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
class object, an `S7_class`.

## Details

Written directly, `S7_inherits(criterion, criterion)` inside a function
whose formal is `criterion` compares the value with itself, and S7 stops
with `` `class` must be an <S7_class> or NULL ``. Reaching the class
through a function of no arguments looks it up in the namespace, past
the formal.
[`check_optimizer_args()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer_args.md)
is the caller.
