# Describe a Line Search to the C++ Side

Flattens a
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
object into the plain list the compiled loop reads, so the C++ side
needs no knowledge of S7. Every search produces the same seven fields
whichever subclass it is, and the compiled code branches on `type`
alone.

## Usage

``` r
line_search_spec(x)
```

## Arguments

- x:

  A
  [`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
  object.

## Value

A list of seven: `type` (`"armijo"` or `"wolfe"`), `c1`, `c2`, `shrink`,
`max_step` (integer), `memory` (integer) and `resolution`.

## Details

A field the given search has no use for is filled with a value the
compiled side ignores: `c2 = 0.9` for the two backtracking searches,
`shrink = 0.5` for Wolfe, `memory = 0` for both of the monotone ones.
One struct reads all three because the shape is fixed.

There are only **two** types.
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
describes itself as `type = "armijo"` with `memory` above zero, the two
differing in the reference value alone, so the compiled loop needs no
third branch. That is the same fact the `memory = 0` identity records
from the other side.
