# Run the Shared Descent Loop

The body of every method that has a direction: prepares the objective,
hands the compiled loop the optimiser's settings and the description of
its direction, and assembles the result.

## Usage

``` r
run_descent(optimizer, fn, par, gr, he, lower, upper, method)
```

## Arguments

- optimizer:

  The
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

- fn, par, gr, he:

  The problem, as the user supplied it.

- lower, upper:

  Box constraints, as in
  [`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

- method:

  A list describing the direction to the compiled loop.

## Value

An
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

## Details

Newton, BFGS, L-BFGS and gradient descent differ only in the `method`
list, which names the direction and carries its parameters. Everything
else — the line search, the stopping rule, the budgets, the trace, the
reporting — is the same code for all of them, which is what makes adding
a fifth method a Direction in C++ and a constructor in R.
