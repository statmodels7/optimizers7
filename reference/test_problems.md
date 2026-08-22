# The Standard Test Problems

The functions optimization papers are argued over: a quadratic, a curved
valley, a few awkward polynomials, two with many minima, and one with a
kink. Each carries its analytic gradient and its known answer.

## Usage

``` r
test_problems(which = NULL)
```

## Arguments

- which:

  An optional character vector naming a subset.

## Value

A named list of problems.

## Details

They are exported for use in testing optimizers generally, not only in
this package's own tests;
[`check_optimizer`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
runs them.

Each element is a list with `name`, `fn`, `gr`, `par` (a starting
point), `solution`, `value`, and two flags. `multimodal` marks a surface
with more than one local minimum, where a local method reaching a
different one is behaving correctly and not failing; `smooth` is `FALSE`
for the one whose derivative does not exist everywhere, where a method
that assumes it does will arrive at the answer and then be unable to
certify it.

The starting points are the ones customarily used, which for
`rosenbrock` and `powell` means the deliberately unhelpful ones the
functions were designed around.

Every problem here has \\f(x^{\ast}) = 0\\, which is a property of the
battery and not of optimization, and it hides one effect. A line search
accepts a step only when the objective decreases by a definite amount,
and near the optimum that decrease is about \\\lVert g \rVert^{2} /
(2\lambda)\\; once it falls below the rounding of the objective itself,
\\\varepsilon \lvert f \rvert\\, no step can be verified and the search
refuses all of them. The smallest gradient a run can reach is therefore
of order

\$\$\lVert g \rVert\_{\text{floor}} \approx \sqrt{2 \lambda \varepsilon
\lvert f(x^{\ast}) \rvert},\$\$

which grows with the value at the solution and is exactly zero for every
problem below. A log-likelihood is the opposite case, being of order one
at its optimum, so an optimizer that reaches `1e-15` here may stop at
`1e-8` there; the defaults of
[`crit_grad`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
allow for that. Adding a constant to any of these objectives moves
neither the minimizer nor the gradient and reproduces the effect.

## See also

[`check_optimizer`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)

## Examples

``` r
names(test_problems())
#> [1] "sphere"     "rosenbrock" "booth"      "beale"      "powell"    
#> [6] "himmelblau" "rastrigin"  "abs_sum"   

p <- test_problems("rosenbrock")[[1]]
minimize(bfgs(), p$fn, p$par, gr = p$gr)@par
#> [1] 0.9999995 0.9999989
p$solution
#> [1] 1 1
```
