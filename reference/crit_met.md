# Has the Stopping Rule Been Met?

Evaluates a criterion on the state of the iteration just completed and
returns whether the run should stop. This is the one generic that a
criterion must implement: a class that inherits from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
and has a method here is a stopping rule, and every algorithm in the
package consults it.

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

`state` has the following components:

- `iter`:

  the iteration just completed.

- `f_new`, `f_old`:

  the objective after and before it.

- `x_new`, `x_old`:

  the parameter vectors, likewise.

- `gradient`:

  the gradient at `x_new`, or `NULL` when the method does not compute
  one.
  [`prox_grad()`](https://statmodels7.github.io/optimizers7/reference/prox_grad.md)
  supplies the proximal gradient mapping instead.

- `stationarity`:

  a non-negative measure of remaining progress, supplied in place of a
  gradient by
  [`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md),
  [`compass()`](https://statmodels7.github.io/optimizers7/reference/compass.md),
  [`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md),
  [`sa()`](https://statmodels7.github.io/optimizers7/reference/sa.md)
  and
  [`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md),
  or `NULL`. See
  [`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md).

The shipped methods supply `f_old` and `x_old` from the first call on,
except
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md),
which supplies them from its second iteration. A descent method calls
the rule with `f_old` and `x_old` set to `NULL` only after a line search
that accepted no step, so that only the state at `x_new` can end the run
there. With box constraints, `x_new`, `x_old` and `gradient` are on the
unconstrained scale on which the optimizer moves.

A rule that needs a component that a method does not supply, such as a
gradient on a derivative-free method, declares it through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
and the optimizer then rejects the rule when the run starts.

## See also

[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md),
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md)

## Examples

``` r
st <- list(iter = 3, f_new = 1.0000001, f_old = 1.0000002,
           x_new = c(1, 2), x_old = c(1, 2), gradient = c(1e-9, -2e-9))
crit_met(crit_grad(1e-8), st)
#> [1] TRUE
crit_met(crit_abs_obj(1e-12), st)
#> [1] FALSE
```
