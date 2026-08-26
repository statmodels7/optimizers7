# How Many Parameters the Objective Takes

Works out the length of the parameter vector by trying lengths and
seeing which the objective accepts. Called by
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
when a starter was given without `npar` and the bounds do not say.

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

What decides it is whether the objective genuinely rejects the wrong
length, and the objective a modeling package hands over usually does:
`X %*% beta` with a parameter of the wrong length is an error rather
than a number, so a regression of any kind is settled at once. A
gradient helps when it spells its components out, since such a gradient
returns a fixed number of them whatever it is handed.

A vectorized objective written in terms of the parameter alone is
another matter: both plausible guesses about R are wrong. Recycling
warns only when the shorter length is not a *divisor* of the longer, so
`sum((p - c(1, 2, 3))^2)` accepts a length-one vector in silence and
returns a perfectly finite 14; and its gradient `2 * (p - c(1, 2, 3))`
returns six components for a length-six argument, so the gradient rule
passes it too. Rosenbrock's `100 (p_2 - p_1^2)^2 + (1 - p_1)^2` accepts
every length from two upwards for the same kind of reason. None of this
is a defect and none of it can be guessed at.

The probe therefore settles the objectives that have a fixed width built
into them and rejects the ones that do not, naming the two lengths it
found. When it rejects, `npar` or a vector of bounds is one word.

The search stops as soon as a *second* length is accepted: the answer is
then already known to be ambiguous and there is no reason to keep
probing. The cost is therefore two evaluations when the objective
accepts any length, and `npar_max` when it accepts exactly one. Either
way it happens once, before the run.

## See also

[`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md)
and
[`start_runif()`](https://statmodels7.github.io/optimizers7/reference/start_runif.md),
the starters that make this question arise, and
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
which asks it.

## Examples

``` r
# A hand-written gradient pins it exactly.
f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                    200 * (p[2] - p[1]^2))
infer_npar(f, gr, function(k) numeric(k))
#> [1] 2

# Without one, the same objective is happy with any length from two
# upwards, and the refusal names the two lengths that decided it.
try(infer_npar(f, NULL, function(k) numeric(k)))
#> Error : The objective accepts more than one length of parameter vector (2 and 3),
#>   so the number of parameters cannot be worked out from it. Say how many,
#>   as in start_zeros(npar = 3).

# An objective with a width built in is settled without a gradient, and
# this is the ordinary case for a model.
set.seed(1)
X <- matrix(rnorm(40 * 3), 40, 3)
y <- as.numeric(X %*% c(1, -2, 0.5) + rnorm(40))
infer_npar(function(b) sum((y - X %*% b)^2), NULL, function(k) numeric(k))
#> [1] 3
```
