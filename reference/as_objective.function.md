# An Ordinary R Function as an Objective

The shipped method: `fn(par)` returns the number to be minimized, and
`gr(par)` and `he(par)` return its gradient and Hessian when the caller
has them. Every evaluation is a callback into R. The compiled objective
counts the evaluations, and the result reports them in `counts`.

## Arguments

- fn:

  A function of the parameter vector returning a single number.

- gr:

  A function of the parameter vector returning the gradient, or `NULL`
  for a central difference of `fn`.

- he:

  A function of the parameter vector returning the Hessian, or `NULL`.
  It is read by
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
  and
  [`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md).

- ...:

  Unused.

## Value

The six-component list described under
[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md),
with `kind = "r"`.

## Validation

`gr` and `he` are each checked to be a function or `NULL`, and any other
value raises an error before the run starts. `fn` is not checked here:
dispatch has already established that it is a function, and the compiled
objective checks at every evaluation that it returns a single number.

## Examples

``` r
obj <- as_objective(function(p) sum(p^2), gr = function(p) 2 * p)
obj$kind
#> [1] "r"
c(obj$has_gradient, obj$has_hessian)
#> [1]  TRUE FALSE

# A gradient that is not a function is rejected before any evaluation.
try(as_objective(function(p) sum(p^2), gr = 1))
#> Error : 'gr' must be a function or NULL.
```
