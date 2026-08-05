# Assemble the Result of a Run

Turns what the C++ loop returned into an
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

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

  The generator state the run began with, or `NULL` for a method that
  draws no random numbers.

## Value

An
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md).

## Details

The one judgement here is the meaning of `converged`: it is taken
straight from whether the stopping rule fired, and never inferred from
the run having ended. An optimizer that exhausted its iterations has not
converged, and saying otherwise turns a failure into a wrong answer that
looks right.
