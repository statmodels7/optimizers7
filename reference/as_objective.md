# Normalize an Objective for the Optimizers

Turns the objective into the single handle the algorithms are written
against.

## Usage

``` r
as_objective(fn, gr = NULL, he = NULL, ...)
```

## Arguments

- fn:

  The objective.

- gr:

  An optional gradient.

- he:

  An optional Hessian. Only Newton uses one; the other methods accept it
  and ignore it, so calling code need not branch on the method.

- ...:

  Passed to methods.

## Value

A list describing the objective to the C++ side: `kind`, the pieces
belonging to that kind, and flags saying which derivatives were supplied
rather than differenced.

## See also

[`crit_met`](https://statmodels7.github.io/optimizers7/reference/crit_met.md),
[`crit_needs`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md),
[`check_criterion`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)

## Examples

``` r
# a plain function, with and without a gradient
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
```
