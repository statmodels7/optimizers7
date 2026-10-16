# Record the Random Number Generator's State

Returns the value of `.Random.seed` before a run that draws random
numbers, so that the run can be repeated exactly.

## Usage

``` r
capture_seed()
```

## Value

An integer vector, the value of `.Random.seed`.

## Details

Assigning the recorded value back into the global environment reproduces
the run. When `.Random.seed` does not exist yet, one uniform is drawn to
make R create it, which costs a single discarded number on the first
stochastic call of a session.

The stream is not created with `set.seed(NULL)`, which reseeds from the
clock and discards the seed that the caller had set.
