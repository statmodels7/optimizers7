# S7 Class for the Gradient Criterion

The rule
[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
builds, and the two methods it implements. It reads `state$gradient`,
takes its max-norm or its 2-norm, and fires when that falls below `tol`.
It declares `"gradient"` through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
so a derivative-free optimizer refuses it when the run starts.

## Usage

``` r
CritGrad(label = character(0), tol = integer(0), norm = character(0))
```

## Arguments

- tol:

  The tolerance, a single positive number.

- norm:

  Either `"max"` or `"2"`.

## Value

An S7 object of class `CritGrad` inheriting from
[`criterion()`](https://statmodels7.github.io/optimizers7/reference/criterion.md),
carrying `label`, `tol` and `norm`.

## Details

[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
returns `FALSE` when the gradient is `NULL`, empty or carries an `NA`,
so a state that has not yet filled it in never accidentally satisfies
the rule. Which norm is used is the `norm` property: on a gradient of
\\(3, 4)\\, `crit_grad(4.5, "max")` fires and `crit_grad(4.5, "2")` does
not.

## See also

[`crit_grad()`](https://statmodels7.github.io/optimizers7/reference/crit_grad.md)
for the constructor and the tolerance it can attain,
[`crit_stationary()`](https://statmodels7.github.io/optimizers7/reference/crit_stationary.md)
for what a derivative-free method reads.

## Examples

``` r
st <- list(iter = 1, f_new = 1, f_old = 2, x_new = 1, x_old = 0,
           gradient = c(3, 4))
crit_needs(crit_grad())
#> [1] "gradient"
c(max = crit_met(crit_grad(4.5, "max"), st),
  two = crit_met(crit_grad(4.5, "2"), st))
#>   max   two 
#>  TRUE FALSE 

# A gradient that is not there is not a small gradient.
crit_met(crit_grad(), list(f_new = 1, gradient = NULL))
#> [1] FALSE
```
