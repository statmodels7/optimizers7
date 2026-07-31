# Turn a Starter Into a Starting Value

The whole of what
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
does with a starter, in one place: settle the number of parameters, draw
the values on the unconstrained scale, and map them back through the
box.

## Usage

``` r
resolve_start(par, fn, gr, lower, upper)
```

## Arguments

- par:

  Whatever was passed as `par`.

- fn, gr:

  The objective and its gradient.

- lower, upper:

  The bounds.

## Value

A numeric vector on the parameter scale.

## Details

The number of parameters is looked for in three places, in order of how
much the caller was willing to say. `npar` on the starter itself is
taken as given. Failing that, a `lower` or `upper` of length greater
than one answers the question, since bounds are one per parameter.
Failing both,
[`infer_npar`](https://statmodels7.github.io/optimizers7/reference/infer_npar.md)
probes the objective.

A numeric `par` passes through untouched, so this costs nothing at all
for the ordinary call.
