# Normalise Box Constraints

Recycles `lower` and `upper` to the length of `par`, checks them, and
returns them in the shape the compiled loop reads.

## Usage

``` r
check_bounds(lower, upper, par)
```

## Arguments

- lower, upper:

  Numeric, of length one or `length(par)`.

- par:

  The starting value.

## Value

A list with one `c(lower, upper)` pair per parameter, or an empty list
when no bound is finite. Call it at the top of a
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
method of your own, then hand each pair to
[`bounded_transform`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md);
an empty list means there is no box and the whole reparametrisation
should be skipped.

## Details

The interface is two numeric vectors rather than a list of pairs, which
is what [`stats::optim`](https://rdrr.io/r/stats/optim.html) and
[`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html) take and what
makes `lower = 0` say "every parameter is positive" in four characters
instead of a list of identical pairs one per coefficient. The list of
pairs is an internal shape, produced here, because that is what the
per-coordinate transform on the C++ side wants.

Recycling is length one or length `p` and nothing between: a `lower` of
length 2 for three parameters is far more likely to be a mistake than a
request, and R's ordinary recycling would silently oblige.

The starting value must be strictly interior, and refusing a boundary
start is not pedantry. The reparametrisation sends a bound to an
infinite value of the transformed variable, so a run started exactly on
one begins at infinity: every subsequent quantity is non-finite and the
failure surfaces far from its cause. Saying so here, naming the
coordinate, costs one check.

## See also

[`bounded_transform`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md),
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

## Examples

``` r
check_bounds(0, Inf, par = c(1, 2))
#> [[1]]
#> [1]   0 Inf
#> 
#> [[2]]
#> [1]   0 Inf
#> 
check_bounds(c(0, -5), c(1, 5), par = c(0.5, 0))
#> [[1]]
#> [1] 0 1
#> 
#> [[2]]
#> [1] -5  5
#> 

# no finite bound is no box at all
length(check_bounds(-Inf, Inf, par = c(1, 2)))
#> [1] 0

# and a start on a bound is refused, naming the coordinate
try(check_bounds(0, 1, par = c(0.5, 1)))
#> Error : The starting value for parameter 2 must lie strictly inside its bounds (0, 1); it is 1.
```
