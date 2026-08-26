# Has the Stopping Rule Been Met?

Asks a criterion whether the run should stop, given the state of the
iteration just completed. This is the one generic a criterion must
implement: a class inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
with a method here is a stopping rule, and every algorithm in the
package will consult it.

## Usage

``` r
crit_met(criterion, state)
```

## Arguments

- criterion:

  A
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object.

- state:

  A named list describing the current iteration; see Details.

## Value

A single logical.

## Details

`state` carries everything any rule could need:

- `iter`:

  the iteration just completed.

- `f_new`, `f_old`:

  the objective after and before it.

- `x_new`, `x_old`:

  the parameter vectors, likewise.

- `gradient`:

  the gradient at `x_new`, or `NULL` when the method does not compute
  one.

- `stationarity`:

  a non-negative measure of remaining progress, supplied by the
  derivative-free methods in place of a gradient, or `NULL`. See
  [`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md).

A rule needing something absent from `state`, such as a gradient on a
derivative-free method, says so through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
so that the optimizer can refuse it when the run starts.

## See also

[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md),
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)

## Examples

``` r
st <- list(iter = 3, f_new = 1.0000001, f_old = 1.0000002,
           x_new = c(1, 2), x_old = c(1, 2), gradient = c(1e-9, -2e-9))
crit_met(crit_grad(1e-8), st)
#> [1] TRUE
crit_met(crit_abs_obj(1e-12), st)
#> [1] FALSE
```
