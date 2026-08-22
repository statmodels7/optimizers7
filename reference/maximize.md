# Maximize a Function

Runs
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
on the negated objective and reports the value with its sign restored.

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
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  object, carrying the algorithm and its settings.

- fn:

  The objective: a function of the parameter vector returning a single
  number.

- par:

  A numeric vector of starting values.

- gr:

  An optional gradient function. Ignored when `fn` carries its own
  gradient.

- he:

  An optional Hessian function. Methods that do not use one ignore it,
  so calling code need not branch on the algorithm.

- lower, upper:

  Box constraints. Numeric, of length one — applying to every parameter
  — or one value per parameter. Default `-Inf` and `Inf`, which is no
  constraint at all. See Details.

- ...:

  Passed to methods.

## Value

An
[`optimizer_result`](https://statmodels7.github.io/optimizers7/reference/optimizer_result.md),
whose `value` and `gradient` refer to the original objective rather than
the negated one.

## Details

\$\$\arg\max\_{x} f(x) = \arg\min\_{x} \\-f(x)\\, \qquad \max\_{x} f(x)
= -\min\_{x}\\-f(x)\\,\$\$

so the point is the one
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
returns and the value, the gradient and the traced objective have their
sign restored.

Every algorithm in the package minimizes, always; that is the convention
and the names say so. This is the thin wrapper for the other direction,
and it exists so that nobody has to remember to negate a log-likelihood
by hand and then negate the answer back.

## See also

[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

## Examples

``` r
maximize(gd(), function(p) -sum((p - c(1, 2))^2), c(0, 0))
#> <optimizer_result> gradient descent
#>   value      : -2.80957e-22
#>   par        : 1 2
#>   iterations : 1   evaluations: f 11, g 0
#>   elapsed    : 1e+03 us
#>   converged  : yes (gradient (max-norm) < 1e-06 or |df| < 1e-10 or |dx| < 1e-08)
#>   note       : gradient obtained by finite differences
```
