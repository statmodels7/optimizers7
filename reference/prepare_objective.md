# Prepare an Objective for an Algorithm

Normalizes the objective, checks that this optimizer can evaluate the
stopping rule, and checks the starting value.

## Usage

``` r
prepare_objective(optimizer, fn, par, gr = NULL, hess = NULL)
```

## Arguments

- optimizer:

  The
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

- fn, gr, hess:

  The objective and its optional derivatives, as supplied.

- par:

  The starting value.

## Value

The objective handle from
[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md).

## Details

A rule that needs a gradient, used with a method that computes none,
would test `NULL` at every iteration and never fire, so the run would
end on its iteration budget and report a failure unrelated to the state
of the run. The check rejects such a rule when the run starts, with an
error that names the missing quantity and the method.
