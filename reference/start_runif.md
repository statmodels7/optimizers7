# Start From a Uniform Draw on the Unconstrained Scale

Each coordinate is drawn independently from `runif(min, max)` on the
unconstrained scale and mapped back through the bounds, so no draw is
ever rejected for being outside the box.

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
  [`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Value

An S7 object of class `UniformStart`.

## Details

The range is in unconstrained units, which is what makes a single
default workable. A draw in \\(-1, 1)\\ becomes a variance between
\\0.37\\ and \\2.7\\, a probability between \\0.27\\ and \\0.73\\, and a
parameter bounded on both sides lands well inside its interval; the same
numbers on the parameter scale would mean quite different things and
would sometimes be inadmissible.

Widen it when the scale of the problem is unknown. `start_runif(-5, 5)`
spans four orders of magnitude for a positive parameter, which is
usually more than enough and is still a range no draw can fall out of.

The draw uses R's ordinary generator, so
[`set.seed`](https://rdrr.io/r/base/Random.html) reproduces it, and the
seed is recorded in the result.

## See also

[`start_zeros`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md),
[`starting_values`](https://statmodels7.github.io/optimizers7/reference/starting_values.md)

## Examples

``` r
f <- function(p) sum((p - c(1, 2, 3))^2)
set.seed(1)
minimize(bfgs(), f, start_runif(npar = 3))@par
#> [1] 1 2 3

# a wider net, and a positive parameter
set.seed(1)
minimize(bfgs(), function(p) (log(p) - 1)^2, start_runif(-5, 5, npar = 1),
         lower = 0)@par
#> [1] 2.718282
```
