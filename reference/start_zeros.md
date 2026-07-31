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
  [`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

## Value

An S7 object of class `ZeroStart`.

## Details

Zero is the right constant only because it is applied on the
unconstrained scale. A vector of zeros on the *parameter* scale is not a
starting point at all for a model with a scale parameter in it: it sits
exactly on the boundary, where the log-likelihood is usually infinite
and the gradient certainly is.

## See also

[`start_runif`](https://statmodels7.github.io/optimizers7/reference/start_runif.md),
[`starting_values`](https://statmodels7.github.io/optimizers7/reference/starting_values.md)

## Examples

``` r
f <- function(p) sum((p - c(1, 2, 3))^2)
minimize(bfgs(), f, start_zeros(3))@par
#> [1] 1 2 3

# zero on the unconstrained scale is one on a positive parameter's scale
minimize(bfgs(), f, start_zeros(3), lower = 0)@par
#> [1] 1 2 3
```
