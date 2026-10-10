# Latin Hypercube Starting Points

Scatters `n - 1` starting points around `par`, generated on the
unconstrained scale so that any bounds are respected by construction.

## Usage

``` r
make_starts(par, n, spread, bounds)
```

## Arguments

- par:

  The starting value of the caller.

- n:

  Total number of starts, including `par`.

- spread:

  Scale of the half-width of the sampling range in unconstrained units:
  the half-width is \\3\\ \texttt{spread}\\ \max(1,
  \lvert\eta_0\rvert)\\.

- bounds:

  Box constraints in the shape
  [`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
  returns, possibly empty.

## Value

A matrix with `n` rows, the first of which is `par`.

## Details

The range of each coordinate is cut into `n - 1` strata, each used
exactly once, so the starts cannot all fall in one corner, as
independent draws can. On the unconstrained scale the range is the
starting value \\\eta_0\\ plus or minus \\3\\ \texttt{spread}\\ \max(1,
\lvert\eta_0\rvert)\\. The unconstrained scale is the parameter itself
for an unbounded coordinate and its shifted log or scaled logit for a
bounded one, so a start for a variance is drawn on the log scale and
comes back positive, without a rejected draw.
