# Translate an Evaluation Budget for the Compiled Loop

The compiled loops hold the budget as an integer, and the default budget
is `Inf`: no cap at all. This maps `Inf` to the largest representable
integer, which no run reaches, and any finite value to itself.

## Usage

``` r
budget_int(x)
```

## Arguments

- x:

  The `max_eval` property, a single positive number or `Inf`.

## Value

A single integer.
