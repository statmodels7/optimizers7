# Check That an Optimiser Keeps Its Promises

Runs an optimiser through a series of checks on what it *reports*, and
then through the standard test problems. Written for whoever adds a
method of their own, and run against every method here.

## Usage

``` r
check_optimizer(
  optimizer,
  problems = test_problems(),
  verbose = TRUE,
  tol = 1e-06
)
```

## Arguments

- optimizer:

  The
  [`optimizer`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  to check.

- problems:

  The battery; defaults to
  [`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md).

- verbose:

  Print the report? Defaults to `TRUE`.

- tol:

  Tolerance for the checks that compare numbers. Defaults to `1e-6`.

## Value

Invisibly, a named list: `checks`, a logical vector with one entry per
numbered check, and `battery`, a data frame of gaps.

## Details

**What this checks is the contract, not the power.** Those are different
questions and conflating them would make the function useless: gradient
descent does not solve Rosenbrock in five hundred iterations, and it is
not broken — it is slow, which is a documented property of the method
and not a defect for a validator to report. So the numbered checks are
all statements an optimiser must satisfy however weak it is, and how
strong it is comes afterwards, as a table of gaps rather than as a
verdict.

The one performance requirement among the numbered checks is that the
optimiser minimises a quadratic. That is a floor no correct method can
fail.

### The checks

1.  `value` is the objective at `par`. A method that reports a value
    from a point it has since left is the kind of defect that survives
    every test written in terms of the value alone.

2.  the reported gradient is the gradient at `par`, checked only for
    optimisers that offer `"gradient"` to a stopping rule and so are
    claiming it is one.
    [`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
    reports an aggregate subgradient and does not make that claim, so it
    is not held to it.

3.  `converged` follows the stopping rule and is never inferred from the
    run having ended. Checked by starving the optimiser of iterations: a
    run cut off after one must not report success.

4.  budgets are respected — `iterations` never exceeds `maxit`.

5.  evaluations are counted; a method reporting zero of them did not
    evaluate anything.

6.  the trace, when kept, is a data frame whose iteration numbers run
    from one and increase.

7.  bounds are respected **strictly**: a probability of exactly 1 is not
    a probability inside \\(0, 1)\\, and the caller's next act is
    usually to divide by it.

8.  the run repeats. A deterministic method must give the same answer
    twice; a stochastic one must give it again from the seed it
    recorded, which tests the recording as well as the repeatability.

9.  [`maximize`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
    is
    [`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
    of the negative.

10. a stopping rule the optimiser cannot evaluate is refused, rather
    than accepted and left never to fire.

11. a starting point where the objective is not finite is an error, not
    a run that quietly returns `NaN`.

12. it minimises a quadratic.

### Reading the battery

The table reports the gap between the value reached and the known
minimum, and it is information rather than judgement. A large gap on
`rastrigin` or `himmelblau` means the method found a different local
minimum, which for a local method is correct behaviour; a large gap on
`abs_sum` means it was defeated by a kink, which is what
[`bundle`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
and the derivative-free methods are for.

## See also

[`test_problems`](https://statmodels7.github.io/optimizers7/reference/test_problems.md),
[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md)

## Examples

``` r
check_optimizer(bfgs())
#> Checking optimizer: BFGS
#>   [ 1] value agrees with par:            [PASSED]
#>   [ 2] gradient agrees with par:         [PASSED]
#>   [ 3] convergence is not assumed:       [PASSED]
#>   [ 4] budgets are respected:            [PASSED]
#>   [ 5] evaluations are counted:          [PASSED]
#>   [ 6] trace is well formed:             [PASSED]
#>   [ 7] bounds are respected strictly:    [PASSED]
#>   [ 8] the run repeats:                  [PASSED]
#>   [ 9] maximize mirrors minimize:        [PASSED]
#>   [10] an unevaluable rule is refused:   [PASSED]
#>   [11] a bad starting point is an error: [PASSED]
#>   [12] it minimises a quadratic:         [PASSED]
#> 
#>   All checks passed.
#> 
#>   battery (gap from the known minimum; information, not a verdict)
#>     sphere       gap  0.00e+00  conv        3 evals  
#>     rosenbrock   gap  1.98e-20  conv       65 evals  
#>     booth        gap  4.80e-18  conv       15 evals  
#>     beale        gap  2.40e-17  conv       21 evals  
#>     powell       gap  3.65e-16  conv       75 evals  
#>     himmelblau   gap  3.04e-19  conv       19 evals  multimodal
#>     rastrigin    gap  7.96e+00  conv       15 evals  multimodal
#>     abs_sum      gap  1.26e-02  -          84 evals  non-smooth

# a method that does not compute a gradient is held to fewer claims, and to
# the same standard on the ones it does make
check_optimizer(nelder_mead(), problems = test_problems("sphere"))
#> Checking optimizer: nelder-mead
#>   [ 1] value agrees with par:            [PASSED]
#>   [ 2] gradient agrees with par:         [PASSED]
#>   [ 3] convergence is not assumed:       [PASSED]
#>   [ 4] budgets are respected:            [PASSED]
#>   [ 5] evaluations are counted:          [PASSED]
#>   [ 6] trace is well formed:             [PASSED]
#>   [ 7] bounds are respected strictly:    [PASSED]
#>   [ 8] the run repeats:                  [PASSED]
#>   [ 9] maximize mirrors minimize:        [PASSED]
#>   [10] an unevaluable rule is refused:   [PASSED]
#>   [11] a bad starting point is an error: [PASSED]
#>   [12] it minimises a quadratic:         [PASSED]
#> 
#>   All checks passed.
#> 
#>   battery (gap from the known minimum; information, not a verdict)
#>     sphere       gap  2.14e-17  conv      276 evals  
```
