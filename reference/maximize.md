# Maximize a Function

Runs
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
on the negated objective and hands back a result whose value, gradient
and traced objective are those of the objective as written. Every
algorithm in the package minimizes; this wrapper maximizes, so that a
log-likelihood need not be negated by hand and the answer negated back.

## Usage

``` r
maximize(
  optimizer,
  fn,
  par,
  gr = NULL,
  hess = NULL,
  lower = -Inf,
  upper = Inf,
  ...
)
```

## Arguments

- optimizer:

  An
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  object carrying the algorithm and its settings.

- fn:

  The objective, a plain function of the parameter vector returning a
  single number, to be maximized. Unlike
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
  an object of another class is rejected and not passed to
  [`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md);
  such an objective has to be negated by the caller and passed to
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

- par:

  A numeric vector of starting values, or a starter object such as
  [`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md).
  With bounds it must lie strictly inside them.

- gr:

  The gradient of `fn` as written, a function of the parameter vector,
  or `NULL` for a central difference. It is negated here, so a caller
  supplies the gradient of the function being maximized.

- hess:

  The Hessian of `fn` as written, or `NULL`. It is negated here too, and
  only
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
  and
  [`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
  read it.

- lower, upper:

  Box constraints, numeric of length 1 or of length `p`, defaulting to
  `-Inf` and `Inf`. Passed through to
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
  unchanged; the sign convention does not touch them.

- ...:

  Passed to
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
  and from there to the method. A `t3` argument for
  [`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
  is the contraction of the third derivative of `fn` as written, and is
  negated here like `gr` and `hess`.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
with the same twelve properties that
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
returns, in which `value`, `gradient` and `trace$value` refer to `fn` as
written.

## Details

\$\$\arg\max\_{x} f(x) = \arg\min\_{x} \\-f(x)\\, \qquad \max\_{x} f(x)
= -\min\_{x}\\-f(x)\\,\$\$

so `par` is exactly what
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
returned, while `value`, `gradient` and the `value` column of `trace`
have their sign restored. The other properties of the result pass
through unchanged: `counts`, `iterations`, `converged`, `criterion_met`
and the optimizer.

The rule named in `criterion_met` was evaluated on the negated
objective. Each shipped rule gives the same verdict on \\f\\ and on
\\-f\\, since each reads a gradient norm, an absolute change or a
relative change.

`fn` here must be a plain function. An objective of some other class,
reached through a method of
[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md),
is rejected with an error: this wrapper negates by composing a closure,
and it cannot negate an object whose evaluation it does not perform.

## See also

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

## Examples

``` r
# A log-likelihood for the normal mean is maximized directly.
set.seed(1)
y <- rnorm(200, mean = 3)
ll <- function(p) sum(dnorm(y, mean = p, log = TRUE))
fit <- maximize(bfgs(), ll, par = 0)
c(fit@par, mean(y))
#> [1] 3.03554 3.03554

# The reported value is the log-likelihood itself, not its negative.
all.equal(fit@value, ll(fit@par))
#> [1] TRUE

# An objective that is not a plain function is rejected.
try(maximize(bfgs(), 1:3, c(0, 0)))
#> Error : maximize() takes a plain function; negate other objectives yourself.
```
