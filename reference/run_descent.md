# Run the Shared Descent Loop

The body of every method that has a direction: it prepares the
objective, passes the settings of the optimizer and the description of
its direction to the compiled loop, and assembles the result.

## Usage

``` r
run_descent(optimizer, fn, par, gr, he, lower, upper, method)
```

## Arguments

- optimizer:

  The
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

- fn, par, gr, he:

  The problem, as the user supplied it.

- lower, upper:

  Box constraints, as in
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

- method:

  A list describing the direction to the compiled loop.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

## Details

[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md),
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md),
[`lbfgs()`](https://statmodels7.github.io/optimizers7/reference/lbfgs.md),
[`gd()`](https://statmodels7.github.io/optimizers7/reference/gd.md),
[`cg()`](https://statmodels7.github.io/optimizers7/reference/cg.md),
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
differ only in the `method` list, which names the direction and carries
its parameters. The line search, the stopping rule, the budgets, the
trace and the reporting are the same code for all of them, so another
method of this kind needs a `Direction` in C++ and a constructor in R.
