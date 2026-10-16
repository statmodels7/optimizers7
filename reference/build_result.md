# Assemble the Result of a Run

Converts the list returned by the compiled loop into an
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md):
it fills in the point and the value, the counts and the trace, sets the
verdict, and composes the message.

## Usage

``` r
build_result(out, optimizer, spec, elapsed, seed = NULL)
```

## Arguments

- out:

  The list from the C++ driver.

- optimizer:

  The optimizer that ran.

- spec:

  The objective handle, consulted for whether the gradient was supplied
  or differenced.

- elapsed:

  Seconds.

- seed:

  The generator state with which the run began, or `NULL` for a method
  that does not draw random numbers.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

## Details

`converged` is the flag returned by the loop, and it is never inferred
from the fact that the run ended. The loop sets it when the stopping
rule fires, when a line search with a declared resolution establishes
that every step fails to decrease the objective by more than that
resolution, and, for
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md),
when the bracket reaches the resolution of the iterate. `criterion_met`
is the label of the rule in the first case and a label of its own in the
others; a run that exhausted a budget reports the budget.
