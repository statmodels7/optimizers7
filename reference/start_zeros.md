# Start From Zero on the Unconstrained Scale

A starting point of all zeros, which after the bound transform is the
sensible middle of every parameter's domain: one for a positive
parameter, one half for a probability, zero for an unbounded one.

## Usage

``` r
start_zeros(npar = NULL)
```

## Arguments

- npar:

  The number of parameters. Defaults to `NULL`, meaning work it out from
  the bounds or from the objective; see
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Value

An S7 object of class
[ZeroStart](https://statmodels7.github.io/optimizers7/reference/ZeroStart-class.md),
inheriting from
[starter](https://statmodels7.github.io/optimizers7/reference/starter-class.md),
to be passed as `par`.

## Details

Zero is the right constant only because it is applied on the
unconstrained scale. What it becomes after the bound transform depends
on the box:

|                       |                    |
|-----------------------|--------------------|
| **bounds**            | **starting value** |
| \\(-\infty, \infty)\\ | 0                  |
| \\(0, \infty)\\       | 1                  |
| \\(0, 1)\\            | 0.5                |
| \\(2, 6)\\            | 4                  |
| \\(-\infty, 5)\\      | 4                  |

A vector of zeros on the *parameter* scale is no starting point at all
for a model with a scale parameter in it: it sits exactly on the
boundary, where the log-likelihood is usually infinite and the gradient
certainly is.
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
refuses such a point for that reason.

## See also

[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md)
for a random start,
[`starting_values()`](https://statmodels7.github.io/optimizers7/reference/starting_values.md)
for the generic that resolves it,
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
for how `npar` is settled.

## Examples

``` r
f <- function(p) sum((p - c(1, 2, 3))^2)
minimize(bfgs(), f, start_zeros(3))@par
#> [1] 1 2 3

# Zero on the unconstrained scale is one on a positive parameter's scale,
# and the midpoint of a two-sided box.
starting_values(start_zeros(), 3)
#> [1] 0 0 0
bounded_transform(c(0, Inf), 0)$h
#> [1] 1
bounded_transform(c(2, 6), 0)$h
#> [1] 4

# The count need not be given when the bounds already say it.
minimize(bfgs(), f, start_zeros(), lower = c(0, 0, 0))@par
#> [1] 1 2 3
```
