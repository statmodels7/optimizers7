# Normalize an Objective for the Optimizers

Packs a function to be minimized, together with its gradient and Hessian
when the caller has them, into the handle on which every algorithm in
the package operates. The result is a named list that contains the
pieces and two flags recording which derivatives were supplied.
`as_objective()` is a generic that dispatches on `fn`, so a caller with
another kind of objective registers one method and every algorithm
accepts that objective.

## Usage

``` r
as_objective(fn, gr = NULL, he = NULL, ...)
```

## Arguments

- fn:

  The objective, a function of the parameter vector returning a single
  number to be minimized. An object of some other class is accepted
  whenever a method has been registered for it; without one, dispatch
  fails and reports the class that it was handed.

- gr:

  The gradient, a function of the parameter vector returning a numeric
  vector of the same length. `NULL`, the default, has the gradient
  differenced from `fn` as above. Any other value is an error.

- he:

  The Hessian, a function of the parameter vector returning a `p x p`
  symmetric matrix. `NULL` is the default.
  [`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
  and
  [`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
  read it; every other method accepts it and does not use it, so calling
  code can pass whatever it has without branching on the method.

- ...:

  Passed to the method that dispatch selects. The shipped method for an
  ordinary function does not read it.

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

## Dispatch of the two generics

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
dispatches on the optimizer, so an algorithm is written once for every
kind of objective. `as_objective()` dispatches on the objective, so each
kind of objective is converted once for every algorithm. Dispatching
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
on both would need one method for each pair of an algorithm and a kind
of objective. The package ships one method, for an ordinary R function.

## A gradient obtained by differencing

The compiled loop differences the objective, once per coordinate,
\$\$\hat{g}\_j = \frac{f(x + h e_j) - f(x - h e_j)}{2h}, \qquad h =
\varepsilon^{1/3} \max(1, \lvert x_j \rvert),\$\$ with \\\varepsilon\\
the machine epsilon, so \\h\\ is about `6.06e-06` for a coordinate of
size one. One gradient then costs \\2p\\ evaluations of the objective,
and its relative error is of the order of \\\varepsilon^{2/3}\\.

A run that differences its gradient records this in the `message` field
of the result.

## A Hessian by differencing the gradient

[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
read a Hessian, and the Hessian is obtained by differencing the gradient
once when the caller does not supply one. With an analytic gradient this
is a single differentiation and costs \\2p\\ gradient evaluations.
Without a gradient, the gradient is itself a difference and the Hessian
is a difference of differences, which costs about \\4p^2\\ evaluations
of the objective and is less accurate. The Hessian is symmetrized before
it is returned, since its two triangles differ by the error of the
differences.

## Points where the objective is not finite

If the objective returns `NA`, `NaN` or an infinity at a point other
than the starting value, the algorithms treat that point as unusable,
reject the trial step and continue, so a line search can try points
outside the region where the objective is defined. A value that is not
finite at the starting value is an error.

## See also

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
for the run itself,
[`newton()`](https://statmodels7.github.io/optimizers7/reference/newton.md)
and
[`chebyshev()`](https://statmodels7.github.io/optimizers7/reference/chebyshev.md)
for the methods that read a Hessian,
[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
for the quantities that a method supplies to a stopping rule.

## Examples

``` r
# The flags record what was supplied, and the compiled loop branches on them.
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

# The example measures the cost of the difference on Rosenbrock from the
# origin. Both runs take
# 21 iterations and reach the same point; the evaluation counts differ.
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

# The message of the run that differenced records it.
a@message
#> [1] "gradient obtained by finite differences"

# A class with no method is rejected, and the error names the class.
try(as_objective(1:3))
#> Error : Can't find method for `as_objective(<integer>)`.
```
