# Describe a Line Search to the C++ Side

Converts a
[`line_search()`](https://statmodels7.github.io/optimizers7/reference/line_search.md)
object into the plain list that the compiled loop reads, so the C++ code
needs no knowledge of S7. Every search produces the same seven fields,
whatever its subclass. The compiled search branches on `type`, and the
descent loop reads `memory` to form the reference value.

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

A field that the search does not use is filled with a value that the
compiled code ignores: `c2 = 0.9` for the two backtracking searches,
`shrink = 0.5` for Wolfe, `memory = 0` for the two monotone ones. The
shape is fixed, so one C++ structure reads all three.

There are only **two** types.
[`nonmonotone()`](https://statmodels7.github.io/optimizers7/reference/nonmonotone.md)
is passed as `type = "armijo"` with its `memory`, since it differs from
[`armijo()`](https://statmodels7.github.io/optimizers7/reference/armijo.md)
in the reference value alone, so the compiled search needs no third
branch. With `memory = 0` the two lists are identical.
