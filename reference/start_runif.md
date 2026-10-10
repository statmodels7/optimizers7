# Start From a Uniform Draw on the Unconstrained Scale

Each coordinate is drawn independently from `runif(min, max)` on the
unconstrained scale and mapped back through the bounds, so every draw
lies inside the box.

## Usage

``` r
start_runif(min = -1, max = 1, npar = NULL)
```

## Arguments

- min, max:

  The range to draw from, in unconstrained units. Both default to a
  width of one either side of zero, and both may be given per parameter
  rather than as a single number.

- npar:

  The number of parameters. Defaults to `NULL`, meaning work it out from
  the bounds or from the objective; see
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Value

An S7 object of class
[UniformStart](https://statmodels7.github.io/optimizers7/reference/UniformStart-class.md),
inheriting from
[starter](https://statmodels7.github.io/optimizers7/reference/starter-class.md),
to be passed as `par`.

## Details

The range is in unconstrained units, which makes a single default
usable. A draw in \\(-1, 1)\\ becomes a variance between `0.368` and
`2.72`, a probability between `0.269` and `0.731`, and a parameter
bounded on both sides lies well inside its interval. The same numbers on
the parameter scale would have different meanings for different
parameters and would sometimes be inadmissible.

A wider range suits a problem of unknown scale. `start_runif(-5, 5)`
spans `0.0067` to `148` for a positive parameter, four orders of
magnitude, and every draw still lies inside the box.

The draw uses the ordinary generator of R, so
[`set.seed()`](https://rdrr.io/r/base/Random.html) reproduces it. The
draw is made before the optimizer starts, so the `seed` that a
stochastic optimizer records does not reproduce the draw, and a
deterministic optimizer records no `seed`.

## See also

[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
for the deterministic start,
[`starting_values()`](https://statmodels7.github.io/optimizers7/reference/starting_values.md)
for the generic that resolves it,
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md),
which uses a starter to generate its own starts.

## Examples

``` r
f <- function(p) sum((p - c(1, 2, 3))^2)
set.seed(1)
minimize(bfgs(), f, start_runif(npar = 3))@par
#> [1] 1 2 3

# What the range means on the parameter scale, for a positive parameter.
range(bounded_transform(c(0, Inf), c(-1, 1))$h)
#> [1] 0.3678794 2.7182818
range(bounded_transform(c(0, Inf), c(-5, 5))$h)
#> [1] 6.737947e-03 1.484132e+02

# A wider net on a positive parameter whose scale is unknown.
set.seed(1)
minimize(bfgs(), function(p) (log(p) - 1)^2, start_runif(-5, 5, npar = 1),
         lower = 0)@par
#> [1] 2.718282

# One range per parameter is allowed; a length that is neither 1 nor npar
# is rejected when the draw is made.
set.seed(1)
starting_values(start_runif(c(-1, -10), c(1, 10)), 2)
#> [1] -0.4689827 -2.5575220
try(starting_values(start_runif(c(-1, -2, -3)), 2))
#> Error : 'min' must have length 1 or 2, one per parameter; it has length 3.
```
