# Normalize an Objective for the Optimizers

Packs a function to be minimized, together with whichever of its
gradient and Hessian the caller has, into the one handle every algorithm
in the package is written against. The result is a named list carrying
the pieces and two flags saying which derivatives were supplied.
`as_objective()` is a generic dispatching on `fn`, so a caller holding
some other kind of objective registers one method and every algorithm
accepts it.

## Usage

``` r
as_objective(fn, gr = NULL, he = NULL, ...)
```

## Arguments

- fn:

  The objective, a function of the parameter vector returning a single
  number to be **minimized**. An object of some other class is accepted
  whenever a method has been registered for it; without one, dispatch
  fails and reports the class it was handed.

- gr:

  The gradient, a function of the parameter vector returning a numeric
  vector of the same length. `NULL`, the default, has the gradient
  differenced from `fn` as above. Anything that is neither a function
  nor `NULL` is an error.

- he:

  The Hessian, a function of the parameter vector returning a `p x p`
  symmetric matrix. `NULL` is the default. Only
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
  reads one; every other method accepts it and ignores it, so calling
  code can pass whatever it has without branching on the method.

- ...:

  Passed to the method dispatched on. The shipped method for an ordinary
  function reads nothing from it.

## Value

A named list of six components describing the objective to the compiled
side:

- `kind`:

  character, `"r"` for an ordinary R function.

- `fn`:

  the objective as supplied.

- `gr`:

  the gradient, or `NULL`.

- `he`:

  the Hessian, or `NULL`.

- `has_gradient`:

  logical, whether `gr` was supplied.

- `has_hessian`:

  logical, whether `he` was supplied.

## Two generics, each dispatching once

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
dispatches on the optimizer, so an algorithm is written once whatever
shape the objective arrived in. `as_objective()` dispatches on the
objective, so the shapes are told apart once whatever algorithm is
running. Dispatching
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
on both would need one method per algorithm per shape.

## What happens when no gradient is supplied

The compiled loop differences the objective, once per coordinate,
\$\$\hat{g}\_j = \frac{f(x + h e_j) - f(x - h e_j)}{2h}, \qquad h =
\varepsilon^{1/3} \max(1, \lvert x_j \rvert),\$\$ with \\\varepsilon\\
the machine epsilon, so \\h\\ is about `6.06e-06` for a coordinate of
size one. One gradient then costs \\2p\\ evaluations of the objective.
The accuracy is what a central difference gives: on the Rosenbrock
function at \\(0.5, 0.5)\\ the differenced gradient is right to a
relative `1.4e-10`.

A run that differences its gradient says so in the result's `message`
field, so a fit is never quietly less exact than it looks.

## The Hessian is the one place two differences compose

Only
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
asks for a Hessian, and when none is supplied it is obtained by
differencing the gradient once. With an analytic gradient that is a
single differentiation and costs \\2p\\ gradient evaluations. With no
gradient either, the gradient is itself a difference and the Hessian is
a difference of differences, which costs about \\4p^2\\ evaluations and
loses accuracy accordingly: one Newton iteration on a quadratic in 20
unknowns takes 1682 evaluations of the objective without a gradient and
42 of the gradient with one. The Hessian is symmetrized before it is
returned, the two triangles differing by the differencing error alone.

## A point the objective refuses

An objective that returns `NA`, `NaN` or an infinity at some point
declares that point unusable. The algorithms treat that as an ordinary
event, back away and continue, so a line search may probe outside the
region where the objective is defined without the run failing.

## See also

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
for the run itself,
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
for the one method that reads a Hessian,
[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
for what a method offers a stopping rule in return.

## Examples

``` r
# The flags record what was supplied, and they are what the loop branches on.
str(as_objective(function(p) sum(p^2)))
#> List of 6
#>  $ kind        : chr "r"
#>  $ fn          :function (p)  
#>  $ gr          : NULL
#>  $ he          : NULL
#>  $ has_gradient: logi FALSE
#>  $ has_hessian : logi FALSE
str(as_objective(function(p) sum(p^2), gr = function(p) 2 * p))
#> List of 6
#>  $ kind        : chr "r"
#>  $ fn          :function (p)  
#>  $ gr          :function (p)  
#>  $ he          : NULL
#>  $ has_gradient: logi TRUE
#>  $ has_hessian : logi FALSE

# What the difference costs, on Rosenbrock from the origin. Both runs take
# 21 iterations and land in the same place; only the bill differs.
f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
                   200 * (p[2] - p[1]^2))

a <- minimize(bfgs(), f, c(0, 0))
b <- minimize(bfgs(), f, c(0, 0), gr = g)
unlist(a@counts)          # 126 objective evaluations, no gradients
#>   f   g   h 
#> 126   0   0 
unlist(b@counts)          # 30 objective evaluations and 24 gradients
#>  f  g  h 
#> 30 24  0 
max(abs(a@par - b@par))   # the two answers agree to 1.5e-08
#> [1] 1.466677e-08

# The run that differenced says so.
a@message
#> [1] "gradient obtained by finite differences"

# A class with no method is refused by name.
try(as_objective(1:3))
#> Error : Can't find method for `as_objective(<integer>)`.
```
