# Has the Stopping Rule Been Met?

The one generic a criterion must implement.

## Usage

``` r
crit_met(criterion, state)
```

## Arguments

- criterion:

  A
  [`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
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
  [`crit_stationary`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md).

A rule that needs something absent from `state` — a gradient, from a
derivative-free method — must say so through
[`crit_needs`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md)
rather than silently never firing.

## See also

[`as_objective`](https://statmodels7.github.io/optimizers7/reference/as_objective.md),
[`crit_needs`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
[`check_criterion`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)

## Examples

``` r
st <- list(iter = 3, f_new = 1.0000001, f_old = 1.0000002,
           x_new = c(1, 2), x_old = c(1, 2), gradient = c(1e-9, -2e-9))
crit_met(crit_grad(1e-8), st)
#> [1] TRUE
crit_met(crit_abs_obj(1e-12), st)
#> [1] FALSE
```
