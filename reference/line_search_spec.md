# Describe a Line Search to the C++ Side

Flattens a
[`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
object into the list the compiled code reads.

## Usage

``` r
line_search_spec(x)
```

## Arguments

- x:

  A
  [`line_search`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object.

## Value

A list with `type`, `c1`, `c2`, `shrink` and `max_step`. Fields a given
search does not use are filled with values the C++ side ignores, so that
the structure is the same shape whichever search it describes.
