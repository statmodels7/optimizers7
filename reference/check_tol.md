# Validate a Tolerance

Checks that the value is a single positive number, so that every
criterion constructor rejects invalid values with the same message. Zero
is rejected as well as a negative value, because a rule with a tolerance
of zero can never fire.

## Usage

``` r
check_tol(tol, nm = "tol")
```

## Arguments

- tol:

  The value supplied.

- nm:

  Its name, for the message. Defaults to `"tol"`.

## Value

Invisibly `TRUE`. Raises an error naming `nm` otherwise.

## Details

It is also called by
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
for `t0`, `t_min`, `t_max` and `qp_tol`, by
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md) for
`alpha0`, `alpha_min`, `alpha_max` and `curv_tol`, by
[`nelder_mead()`](https://statmodels7.github.io/optimizers7/reference/nelder_mead.md)
for `degenerate_tol` and by
[`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
for `spread` and `distinct_tol`. These callers pass the name of their
own argument in `nm`, so that the message names the argument that was
supplied.
