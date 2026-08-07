# Stop When the Gradient Is Small

The rule \\\lVert \nabla f \rVert \< \texttt{tol}\\.

## Usage

``` r
crit_grad(tol = 1e-06, norm = c("max", "2"))
```

## Arguments

- tol:

  Numeric tolerance. Defaults to `1e-6`.

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
`1e-6` in the 2-norm is a far stricter demand in high dimension. Both
are available, so the choice is only a default.

How small a gradient a run can actually reach is set by the objective,
not by the method. A line search accepts a step only when the objective
decreases by a definite amount, and near a minimum that decrease is
about \\\lVert \nabla f \rVert^{2} / (2\lambda)\\ for a curvature
\\\lambda\\. Once it drops below the rounding of the objective itself,
about \\\varepsilon \lvert f \rvert\\, no step in any direction can be
verified and the search stops, so the smallest attainable gradient is
around \\\sqrt{2 \lambda \varepsilon \lvert f^{\*} \rvert}\\ and grows
with the value at the solution. On conjugate gradients applied to
Rosenbrock, adding a constant to the objective — which moves neither the
minimizer nor the gradient — takes the attainable gradient from `1.9e-9`
at \\f^{\*} = 0\\ to `4.4e-8` at \\f^{\*} = 1\\ and `6.5e-5` at \\f^{\*}
= 10^{6}\\. The default suits an objective of order one at its solution,
which is what a log-likelihood per observation is; an objective that
lands in the millions needs a correspondingly looser tolerance, and one
that lands at zero can be asked for much more.

Only usable by a method that computes a gradient; a derivative-free
optimizer rejects it rather than accepting a rule that can never fire.

## See also

[`crit_rel_obj`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md),
[`crit_any`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)

## Examples

``` r
crit_grad()
#> <criterion> gradient (max-norm) < 1e-06
crit_grad(1e-10, norm = "2")
#> <criterion> gradient (2-norm) < 1e-10
```
