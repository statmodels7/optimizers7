# Print the Report of check_optimizer

Writes the twelve verdicts, one per line as `[PASSED]` or `[FAILED]`,
then a summary naming every failing check, then the battery as one line
per problem with its gap, its convergence flag, its evaluation count and
its note.

## Usage

``` r
print_optimizer_check(optimizer, ok, battery)
```

## Arguments

- optimizer:

  The
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  checked, read for its name.

- ok:

  The named logical vector of twelve checks.

- battery:

  The data frame
  [`run_battery()`](https://statmodels7.github.io/optimizers7/reference/run_battery.md)
  returned.

## Value

Invisibly `NULL`. Called for the output.
