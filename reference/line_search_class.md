# The line_search Class Object

Returns the class object from the namespace at the time of the call, so
that a class check does not depend on whether the class object is the
one created when the package was loaded.

## Usage

``` r
line_search_class()
```

## Value

The
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
class object.

## Details

Comparing S7 classes is a comparison of objects, so it is `FALSE` for a
class rebuilt from the same definition, which happens under any loader
that re-evaluates the code instead of loading it, covr among them.
