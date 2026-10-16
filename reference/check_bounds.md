# Normalize Box Constraints

Recycles `lower` and `upper` to the length of `par`, checks them, and
returns them in the shape that the compiled loop reads.

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
when every bound is infinite. Call it at the top of a user-defined
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
method, then hand each pair to
[`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md);
an empty list means there is no box and the whole reparametrization
should be skipped.

## Details

The interface is two numeric vectors and not a list of pairs, as in
[`stats::optim`](https://rdrr.io/r/stats/optim.html) and
[`stats::nlminb`](https://rdrr.io/r/stats/nlminb.html), so that
`lower = 0` makes every parameter positive without one pair per
coefficient. The list of pairs is an internal shape, produced here,
because the transform of each coordinate in the compiled code reads it.

Recycling is to length one or length `p` and nothing between: a `lower`
of length 2 for three parameters is much more likely to be a mistake
than a request, and the ordinary recycling of R would accept it
silently.

The starting value must be strictly interior. The reparametrization
sends a bound to an infinite value of the transformed variable, so a run
that starts exactly on a bound would begin at infinity, every later
quantity would be non-finite, and the failure would appear far from its
cause. The error raised here names the coordinate.

## See also

[`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md),
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

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

# infinite bounds on every parameter give no box at all
length(check_bounds(-Inf, Inf, par = c(1, 2)))
#> [1] 0

# and a start on a bound is rejected, naming the coordinate
try(check_bounds(0, 1, par = c(0.5, 1)))
#> Error : The starting value for parameter 2 must lie strictly inside its bounds (0, 1); it is 1.
```
