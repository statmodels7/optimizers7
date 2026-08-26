# Run an Optimizer Over the Battery

Runs the optimizer on each problem from that problem's own starting
point and records the gap between the value reached and the known
minimum. The result is information: no row is a pass or a failure.

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

  A list in the shape
  [`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
  returns.

## Value

A data frame with one row per problem and the columns `problem`
(character), `value` and `gap` (numeric), `converged` (logical),
`evaluations` (integer) and `note` (character, `"multimodal"`,
`"non-smooth"`, an error message, or empty).

## Details

A problem the optimizer cannot run at all is caught rather than
propagated, so one method that refuses one problem does not lose the
other seven. Such a row carries `NA` in every numeric column and the
error message in `note`.
