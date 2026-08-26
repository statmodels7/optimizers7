# Maximize a Function

Runs
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
on the negated objective and hands back a result whose value, gradient
and traced objective are those of the objective as written. Every
algorithm in the package minimizes; this is the wrapper for the other
direction, so that a log-likelihood need not be negated by hand and the
answer negated back.

## Usage

``` r
maximize(
  optimizer,
  fn,
  par,
  gr = NULL,
  he = NULL,
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

  The objective, a **plain function** of the parameter vector returning
  a single number, to be maximized. Unlike
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
  an object of another class is refused rather than passed to
  [`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md);
  negate such an objective yourself and call
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).

- par:

  A numeric vector of starting values, or a starter object such as
  [`start_zeros()`](https://statmodels7.github.io/optimizers7/reference/start_zeros.md).
  With bounds it must lie strictly inside them.

- gr:

  The gradient of `fn` as written, a function of the parameter vector,
  or `NULL` for a central difference. It is negated here, so a caller
  supplies the gradient of the function being maximized.

- he:

  The Hessian of `fn` as written, or `NULL`. Negated here too, and read
  only by
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md).

- lower, upper:

  Box constraints, numeric of length 1 or of length `p`, defaulting to
  `-Inf` and `Inf`. Passed through to
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
  unchanged; the sign convention does not touch them.

- ...:

  Passed to
  [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
  and from there to the method.

## Value

An
[`optimizer_result()`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md)
with the same twelve properties
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
returns, and with `value`, `gradient` and `trace$value` referring to
`fn` as written.

## Details

\$\$\arg\max\_{x} f(x) = \arg\min\_{x} \\-f(x)\\, \qquad \max\_{x} f(x)
= -\min\_{x}\\-f(x)\\,\$\$

so `par` is exactly what
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
returned, while `value`, `gradient` and the `value` column of `trace`
have their sign restored. Everything else on the result belongs to the
run and passes through untouched: `counts`, `iterations`, `converged`,
`criterion_met` and the optimizer.

A consequence worth knowing: `criterion_met` names a rule that was
evaluated on the **negated** objective. For a rule reading a gradient
norm or an absolute change that makes no difference, both being
invariant under the sign. For a rule reading the objective's own scale,
such as
[`crit_rel_obj()`](https://statmodels7.github.io/optimizers7/reference/crit_rel_obj.md),
it makes none either, the relative change being a ratio.

`fn` here must be a plain function. An objective of some other class,
reached through a method of
[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md),
is refused by name: this wrapper negates by composing a closure, and it
has no way to negate an object whose evaluation it does not perform.

## See also

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

## Examples

``` r
# A log-likelihood: the normal mean, maximized directly.
set.seed(1)
y <- rnorm(200, mean = 3)
ll <- function(p) sum(dnorm(y, mean = p, log = TRUE))
fit <- maximize(bfgs(), ll, par = 0)
c(fit@par, mean(y))
#> [1] 3.03554 3.03554

# The reported value is the log-likelihood itself, not its negative.
all.equal(fit@value, ll(fit@par))
#> [1] TRUE

# An objective that is not a plain function is refused by name.
try(maximize(bfgs(), 1:3, c(0, 0)))
#> Error : maximize() takes a plain function; negate other objectives yourself.
```
