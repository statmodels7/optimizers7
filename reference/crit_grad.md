# Stop When the Gradient Is Small

Builds the rule \\\lVert \nabla f \rVert \< \texttt{tol}\\, the test for
a stationary point and the one a maximum likelihood fit should ask for.
Usable only by a method that computes a gradient; a derivative-free
optimizer refuses it when the run starts instead of accepting a rule
that can never fire.

## Usage

``` r
crit_grad(tol = 1e-06, norm = c("max", "2"))
```

## Arguments

- tol:

  Numeric tolerance, a single positive number. Defaults to `1e-6`; see
  below for why that rather than something smaller.

- norm:

  `"max"` (default) or `"2"`. Partial matching applies and any other
  string is refused, naming both.

## Value

An S7 object of class
[CritGrad](https://statmodels7.github.io/optimizers7/reference/CritGrad-class.md),
inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md).

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
Rosenbrock, adding a constant to the objective (which moves neither the
minimizer nor the gradient) takes the attainable gradient from `1.9e-9`
at \\f^{\*} = 0\\ to `4.4e-8` at \\f^{\*} = 1\\ and `6.5e-5` at \\f^{\*}
= 10^{6}\\. The default suits an objective of order one at its solution,
as a log-likelihood per observation is; an objective that lands in the
millions needs a correspondingly looser tolerance, and one that lands at
zero can be asked for much more.

Only usable by a method that computes a gradient; a derivative-free
optimizer rejects it rather than accepting a rule that can never fire.

## See also

[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for what a derivative-free method reads,
[`crit_any()`](https://statmodels7.github.io/optimizers7/reference/crit_any.md)
for the disjunction the gradient methods default to,
[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
for the battery the attainable-gradient figures come from.

## Examples

``` r
crit_grad()
#> <criterion> gradient (max-norm) < 1e-06
crit_grad(1e-10, norm = "2")
#> <criterion> gradient (2-norm) < 1e-10

# The two norms differ, and the max-norm is the looser of the two.
st <- list(f_new = 1, f_old = 2, x_new = 1, x_old = 0, gradient = c(3, 4))
c(max = crit_met(crit_grad(4.5, "max"), st),
  two = crit_met(crit_grad(4.5, "2"), st))
#>   max   two 
#>  TRUE FALSE 

# Asking for the rule alone, where the shipped default is a disjunction.
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))
r <- minimize(bfgs(criterion = crit_grad()), f, c(-1.2, 1), gr = g)
c(r@criterion_met, max(abs(g(r@par))))
#> [1] "gradient (max-norm) < 1e-06" "3.24050564146347e-10"       

# And the constant added to the objective, which moves neither the
# minimizer nor the gradient, changes what is reachable.
shifted <- function(p) f(p) + 1e6
s <- minimize(cg(criterion = crit_grad(1e-14), maxit = 20000), shifted,
              c(-1.2, 1), gr = g)
c(converged = s@converged, gradient = max(abs(g(s@par))))
#>    converged     gradient 
#> 0.000000e+00 6.458007e-05 
```
