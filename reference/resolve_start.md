# Turn a Starter Into a Starting Value

Everything that
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
does with a starter, in one place: it settles the number of parameters,
draws the values on the unconstrained scale, and maps them back through
the box.

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

The number of parameters is looked for in three places, in order. The
`npar` of the starter is taken as given. Otherwise a `lower` or `upper`
of length greater than one determines it, since bounds are one per
parameter. Otherwise
[`infer_npar()`](https://statmodels7.github.io/optimizers7/reference/infer_npar.md)
probes the objective.

A numeric `par` passes through unchanged, so the ordinary call has no
added cost.
