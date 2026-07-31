# Stop When the Gradient Is Small

The rule \\\lVert \nabla f \rVert \< \texttt{tol}\\.

## Usage

``` r
crit_grad(tol = 1e-08, norm = c("max", "2"))
```

## Arguments

- tol:

  Numeric tolerance. Defaults to `1e-8`.

- norm:

  `"max"` (default) or `"2"`.

## Value

A
[`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
object.

## Details

The max-norm is the default because it does not grow with the dimension
the way the 2-norm does: the same tolerance then means the same thing
for a two-parameter problem and a two-hundred-parameter one, whereas
`1e-8` in the 2-norm is a far stricter demand in high dimension. Both
are available, so the choice is only a default.

Only usable by a method that computes a gradient; a derivative-free
optimiser refuses it rather than accepting a rule that can never fire.

## See also

[`crit_rel_obj`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md),
[`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)

## Examples

``` r
crit_grad()
#> <criterion> gradient (max-norm) < 1e-08
crit_grad(1e-10, norm = "2")
#> <criterion> gradient (2-norm) < 1e-10
```
