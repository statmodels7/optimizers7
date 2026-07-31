# Latin Hypercube Starting Points

Scatters `n - 1` starting points around `par`, generated on the
unconstrained scale so that any bounds are respected by construction.

## Usage

``` r
make_starts(par, n, spread, bounds)
```

## Arguments

- par:

  The user's starting value.

- n:

  Total number of starts, including `par`.

- spread:

  Half-width of the sampling range, in unconstrained units.

- bounds:

  Box constraints in the shape
  [`check_bounds`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
  returns, possibly empty.

## Value

A matrix with `n` rows, the first of which is `par`.

## Details

Each coordinate's range is cut into `n - 1` strata, each used exactly
once, so the starts cannot all fall in one corner the way independent
draws can. The range is `par` plus or minus `3 * spread` on the
unconstrained scale, which for an unbounded parameter is the parameter
itself and for a bounded one is its log or logit — so a start for a
variance is drawn as a log and comes back positive without a single
rejected draw.
