# Validate the Settings Every Optimizer Shares

The checks each constructor would otherwise repeat, in one place and in
one wording.

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

Invisibly `TRUE`; raises an error otherwise.
