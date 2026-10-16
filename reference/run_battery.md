# Run an Optimizer Over the Battery

Runs the optimizer on each problem from that problem's own starting
point and records the gap between the value reached and the known
minimum. The result is information: a row is neither a pass nor a
failure.

## Usage

``` r
run_battery(optimizer, problems)
```

## Arguments

- optimizer:

  The
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  to run.

- problems:

  A list in the shape that
  [`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
  returns.

## Value

A data frame with one row per problem and the columns `problem`
(character), `value` and `gap` (numeric), `converged` (logical),
`evaluations` (integer) and `note` (character, `"multimodal"`,
`"non-smooth"`, an error message, or empty).

## Details

An error on a problem is caught and reported in the `note` column, so
that it does not discard the results on the other problems. Such a row
keeps the name of the problem and has `NA` in the other numeric and
logical columns.
