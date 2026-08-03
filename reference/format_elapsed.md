# Format a Duration With a Unit Matched to Its Size

Renders a time in seconds using the unit its magnitude calls for:
microseconds below a millisecond, milliseconds below a second, seconds
below a minute, minutes and seconds below an hour, hours and minutes
above.

## Usage

``` r
format_elapsed(sec)
```

## Arguments

- sec:

  A single non-negative number of seconds.

## Value

A character string, or `NA_character_` when `sec` is missing or not
finite.
