# Validate the Settings Every Optimizer Shares

Checks the seven arguments common to every constructor, so that all
fifteen reject the same nonsense in the same words. Called at the top of
each constructor, before any algorithm-specific validation.

## Usage

``` r
check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
```

## Arguments

- criterion:

  The stopping rule.

- maxit, max_eval, refresh:

  Numeric budgets.

- verbose, keep_trace:

  Logical flags.

## Value

Invisibly `TRUE`. Raises an error naming the offending argument
otherwise.

## Details

The rules, and they differ from one another:

- `criterion` must inherit from
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

- `maxit` must be a single number at least 1 **and finite**.

- `max_eval` must be a single number at least 1; `Inf` passes.

- `refresh` must be a single number at least 0.

- `verbose` and `keep_trace` must each be `TRUE` or `FALSE`, `NA`
  refused.
