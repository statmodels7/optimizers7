# The criterion Class Object

Fetched rather than captured, so that the check above cannot be fooled
by the class being re-created.

## Usage

``` r
criterion_class()
```

## Value

The
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
class object.

## Details

Comparing S7 classes with
[`identical()`](https://rdrr.io/r/base/identical.html) is object
identity, and a class rebuilt from the same definition is not identical
to the original. Under covr, which re-evaluates the code instead of
loading it, that turned every numerical fallback into the chain of first
differences the design exists to forbid – and the local suite,
`R CMD check --as-cran` and a five-platform matrix all passed. Only the
coverage job failed.
