# Record the Random Number Generator's State

The state `.Random.seed` held before a run that draws random numbers, so
that the run can be repeated exactly.

## Usage

``` r
capture_seed()
```

## Value

An integer vector, the value of `.Random.seed`.

## Details

Assigning the recorded value back into the global environment reproduces
the run. When no stream exists yet one uniform is drawn to force R to
create one, which costs a single discarded number on the first
stochastic call of a session.

It is emphatically **not** done with `set.seed(NULL)`, which reseeds
from the clock and throws away whatever the caller had set. That is a
defect this toolkit has already met once: it made a check in
distributions7 silently random, and the resulting test failed about one
CI run in several hundred, on whichever platform happened to draw it.
