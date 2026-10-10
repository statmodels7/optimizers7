# Format a Duration With a Unit Matched to Its Size

Renders a time in seconds in a unit matched to its size: microseconds
below a millisecond, milliseconds below a second and seconds below a
minute, each to three significant figures, then whole minutes and
seconds below an hour, and whole hours and minutes above. The unit is
chosen on the rounded value, so that no field reaches its carry
(`"59 min 60 s"`). Used by
[`print.optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/print.optimizer_result.md)
for the `elapsed` line.

## Usage

``` r
format_elapsed(sec)
```

## Arguments

- sec:

  A single number of seconds.

## Value

A character string such as `"250 ms"`, `"1 min 30 s"` or `"1 h 7 min"`.
`NA_character_` when `sec` is empty, missing or not finite, so an
unmeasured duration is reported as unmeasured.
