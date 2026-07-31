# The Standard Test Problems

The functions optimisation papers are argued over: a quadratic, a curved
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

They are exported because they are useful to anyone writing an
optimiser, not only to this package's own tests:
[`check_optimizer`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
runs them, and so can you.

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

## See also

[`check_optimizer`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)

## Examples

``` r
names(test_problems())
#> [1] "sphere"     "rosenbrock" "booth"      "beale"      "powell"    
#> [6] "himmelblau" "rastrigin"  "abs_sum"   

p <- test_problems("rosenbrock")[[1]]
minimize(bfgs(), p$fn, p$par, gr = p$gr)@par
#> [1] 1 1
p$solution
#> [1] 1 1
```
