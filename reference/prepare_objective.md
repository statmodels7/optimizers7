# Prepare an Objective for an Algorithm

Normalises the objective, checks the stopping rule can actually be
evaluated by this optimiser, and checks the starting value.

## Usage

``` r
prepare_objective(optimizer, fn, par, gr = NULL, he = NULL)
```

## Arguments

- optimizer:

  The
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md).

- fn, gr, he:

  The objective and its optional derivatives, as supplied.

- par:

  The starting value.

## Value

The objective handle from
[`as_objective`](https://statmodels7.github.io/optimizers7/reference/as_objective.md).

## Details

The criterion check is the interesting one. A rule needing a gradient
handed to a method that computes none would sit there testing `NULL` at
every iteration and never fire, so the run would end on the iteration
budget and report failure for a reason nowhere near the truth. Refusing
it here, by name, is the same discipline as `check_link()` in
linkfunctions7 reporting a numerical derivative order as numerical
rather than as passed: a check that cannot be evaluated must say so
rather than pass or fail silently.
