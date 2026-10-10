# How Many Parameters the Objective Takes

Determines the length of the parameter vector by trying lengths and
recording which ones the objective accepts. Called by
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
when a starter was given without `npar` and the bounds do not determine
it.

## Usage

``` r
infer_npar(fn, gr, probe, npar_max = 50)
```

## Arguments

- fn:

  The objective.

- gr:

  Its gradient, or `NULL`. Supplying one makes the answer much more
  likely to be unique; see Details.

- probe:

  A function of one integer returning a candidate parameter vector of
  that length, so that the objective is probed where it will be used.

- npar_max:

  The largest length tried. Defaults to `50`.

## Value

A single integer, the one length accepted. Raises an error naming the
two lengths it found when more than one is accepted, and an error when
none is.

## Details

A length is accepted when `fn` returns a single finite number for it and
raises neither an error nor a warning, and, if `gr` was supplied, when
the gradient comes back with the same length as its argument.

The probe works when the objective rejects a wrong length, as the
objective of a modeling package usually does: `X %*% beta` with a
parameter of the wrong length is an error and not a number, so a
regression model is settled at once. A gradient helps when it writes its
components out, since such a gradient returns a fixed number of them,
whatever its argument.

A vectorized objective written in terms of the parameter alone behaves
differently. Recycling warns only when the shorter length is not a
*divisor* of the longer, so `sum((p - c(1, 2, 3))^2)` accepts a vector
of length one without a warning and returns the finite value 14; and its
gradient `2 * (p - c(1, 2, 3))` returns six components for an argument
of length six, so the gradient test accepts that length too. The
Rosenbrock function `100 (p_2 - p_1^2)^2 + (1 - p_1)^2`, written with
`p[1]` and `p[2]`, accepts every length from two upwards for a related
reason. The length cannot be determined from such an objective.

The probe therefore settles the objectives that have a fixed width built
into them and rejects the others, with an error that names two lengths
that were accepted. The caller then gives the count through `npar` or
through a vector of bounds.

The search stops as soon as a *second* length is accepted, since the
answer is then known to be ambiguous. The cost is therefore two
evaluations when the objective accepts any length, and `npar_max` when
it accepts exactly one. Either way it happens once, before the run.

## See also

[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
and
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md),
the starters that make this determination necessary, and
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
which calls it.

## Examples

``` r
# A hand-written gradient pins it exactly.
f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))
infer_npar(f, gr, function(k) numeric(k))
#> [1] 2

# Without one, the same objective accepts any length from two upwards,
# and the error names the two lengths that were accepted first.
try(infer_npar(f, NULL, function(k) numeric(k)))
#> Error : The objective accepts more than one length of parameter vector (2 and 3),
#>   so the number of parameters cannot be worked out from it. Say how many,
#>   as in start_zeros(npar = 3).

# An objective with a fixed width is settled without a gradient, which is
# the ordinary case for a model.
set.seed(1)
X <- matrix(rnorm(40 * 3), 40, 3)
y <- as.numeric(X %*% c(1, -2, 0.5) + rnorm(40))
infer_npar(function(b) sum((y - X %*% b)^2), NULL, function(k) numeric(k))
#> [1] 3
```
