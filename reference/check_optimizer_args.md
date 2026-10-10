# Validate the Settings Every Optimizer Shares

Checks the six arguments common to the constructors (the stopping rule,
`maxit`, `max_eval`, `verbose`, `refresh` and `keep_trace`), so that
every constructor that calls it rejects invalid settings with the same
messages. It is called at the top of each such constructor, before any
validation specific to the algorithm.

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

The rules are:

- `criterion` must inherit from
  [`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

- `maxit` must be a single number at least 1 **and finite**.

- `max_eval` must be a single number at least 1; `Inf` passes.

- `refresh` must be a single number at least 0.

- `verbose` and `keep_trace` must each be `TRUE` or `FALSE`; `NA` is
  rejected.
