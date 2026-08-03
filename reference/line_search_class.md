# The line_search Class Object

Fetched rather than captured, so that a check cannot be fooled by the
class being re-created.

## Usage

``` r
line_search_class()
```

## Value

The
[`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
class object.

## Details

Comparing S7 classes by identity is object identity, so it is `FALSE`
for a class rebuilt from the same definition – which is what happens
under any loader that re-evaluates the code rather than loading it, covr
among them. The same defect in linkfunctions7 silently turned every
numerical fallback into a chain of first differences, and only the
coverage job noticed.
