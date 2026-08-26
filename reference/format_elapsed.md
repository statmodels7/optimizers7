# Format a Duration With a Unit Matched to Its Size

Renders a time in seconds using the unit its magnitude calls for, to
three significant figures: microseconds below a millisecond,
milliseconds below a second, seconds below a minute, whole minutes and
seconds below an hour, whole hours and minutes above. Used by
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
