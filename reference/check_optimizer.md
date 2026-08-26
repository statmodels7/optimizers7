# Check That an Optimizer Keeps Its Promises

Puts an optimizer through twelve checks on what it *reports*, prints the
verdict on each, and then runs it over the standard test problems and
reports the gap it reached on every one. Written for whoever adds a
method of their own; all twelve shipped methods pass all twelve checks.

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
  [`optimizer()`](https://statmodels7.github.io/optimizers7/reference/optimizer.md)
  to check. Anything else raises an error naming
  [`bfgs()`](https://statmodels7.github.io/optimizers7/reference/bfgs.md)
  as an example.

- problems:

  A list of problems in the shape
  [`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
  returns. Defaults to all eight. Pass a subset to keep the report
  short: `test_problems("sphere")`.

- verbose:

  `TRUE` (default) prints the twelve verdicts and the battery table;
  `FALSE` returns the same information silently.

- tol:

  Tolerance for the checks that compare numbers, a single positive
  number, default `1e-6`. It governs checks 1, 2 and 9; check 12 uses a
  fixed `1e-3` on the distance to the solution, and check 7 tests the
  bounds strictly and reads no tolerance at all.

## Value

Invisibly, a named list of two:

- `checks`:

  a named logical vector of length 12, one entry per numbered check,
  named by what the check asserts.

- `battery`:

  a data frame with one row per problem and the columns `problem`,
  `value`, `gap`, `converged`, `evaluations` and `note`. A problem the
  optimizer could not run at all has `NA` throughout and the error
  message in `note`.

## The contract and the power are different questions

The twelve numbered checks are statements an optimizer must satisfy
however weak it is. How strong it is comes afterwards, as a table of
gaps with no verdict attached. Conflating the two would make the
function useless: gradient descent does not solve Rosenbrock in five
hundred iterations, and it is slow rather than broken, which is a
documented property of the method.

The one performance requirement among the numbered checks is check 12,
minimizing a quadratic. That is a floor no correct method can fail.

## The twelve checks

1.  `value` is the objective at `par`. A method reporting a value from a
    point it has since left is a defect that survives every test written
    in terms of the value alone.

2.  the reported gradient is the gradient at `par`. Checked only for an
    optimizer that offers `"gradient"` through
    [`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md),
    and so is claiming that what it reports is one.
    [`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
    reports an aggregate subgradient and makes no such claim, so it is
    exempt.

3.  `converged` follows the stopping rule and is never inferred from the
    run having ended. Checked by starving the optimizer of iterations: a
    run cut off after one must not report success.

4.  budgets are respected: `iterations` never exceeds `maxit`.

5.  evaluations are counted. A method reporting zero of them evaluated
    nothing.

6.  the trace, when kept, is a data frame whose iteration numbers start
    at one and increase.

7.  bounds are respected **strictly**. A probability of exactly 1 is not
    a probability inside \\(0, 1)\\, and the caller's next act is
    usually to divide by it.

8.  the run repeats. A deterministic method gives the same answer twice;
    a stochastic one gives it again from the seed it recorded, which
    tests the recording as well as the repeatability.

9.  [`maximize()`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
    is
    [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
    of the negative.

10. a stopping rule the optimizer cannot evaluate is rejected at
    construction, not accepted and left never to fire.

11. a starting point where the objective is not finite raises an error,
    so no run quietly returns `NaN`.

12. it minimizes a quadratic.

A deliberately lying optimizer, one returning a made-up point with
`converged = TRUE`, fails six of the twelve: 1, 2, 3, 10, 11 and 12. The
six it passes are the bookkeeping ones, which is why the battery is not
the whole check.

## The problem battery

The table reports the gap between the value reached and the known
minimum, as information. A large gap on `rastrigin` or `himmelblau`
means the method found a different local minimum, which for a local
method is correct behavior, and the `note` column labels those two. A
large gap on `abs_sum` means the method was defeated by a kink, the case
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
and the derivative-free methods exist for.

## See also

[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
for the battery,
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
and
[`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
for the two pieces a method of your own has to use,
[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
for the claim check 2 rests on.

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
#>   [10] an unevaluable rule is rejected:  [PASSED]
#>   [11] a bad starting point is an error: [PASSED]
#>   [12] it minimizes a quadratic:         [PASSED]
#> 
#>   All checks passed.
#> 
#>   battery (gap from the known minimum; information, not a verdict)
#>     sphere       gap  4.62e-33  conv        3 evals  
#>     rosenbrock   gap  3.23e-13  conv       49 evals  
#>     booth        gap  1.73e-17  conv        8 evals  
#>     beale        gap  3.18e-15  conv       17 evals  
#>     powell       gap  4.28e-11  conv       39 evals  
#>     himmelblau   gap  6.27e-14  conv       17 evals  multimodal
#>     rastrigin    gap  0.00e+00  conv       11 evals  multimodal
#>     abs_sum      gap  1.26e-02  -          84 evals  non-smooth

# A method that computes no gradient is held to fewer claims, and to the
# same standard on the ones it does make.
res <- check_optimizer(nelder_mead(), problems = test_problems("sphere"))
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
#>   [10] an unevaluable rule is rejected:  [PASSED]
#>   [11] a bad starting point is an error: [PASSED]
#>   [12] it minimizes a quadratic:         [PASSED]
#> 
#>   All checks passed.
#> 
#>   battery (gap from the known minimum; information, not a verdict)
#>     sphere       gap  2.14e-17  conv      276 evals  
all(res$checks)
#> [1] TRUE

# The battery is a table, so it can be read rather than printed.
check_optimizer(cg(), verbose = FALSE)$battery
#>      problem        value          gap converged evaluations       note
#> 1     sphere 0.000000e+00 0.000000e+00      TRUE           3           
#> 2 rosenbrock 7.399642e-10 7.399642e-10      TRUE         249           
#> 3      booth 1.010522e-11 1.010522e-11      TRUE          73           
#> 4      beale 1.116830e-11 1.116830e-11      TRUE         100           
#> 5     powell 1.048119e-07 1.048119e-07      TRUE         549           
#> 6 himmelblau 1.754022e-11 1.754022e-11      TRUE          97 multimodal
#> 7  rastrigin 0.000000e+00 0.000000e+00      TRUE          60 multimodal
#> 8    abs_sum 4.504423e-02 4.504423e-02      TRUE         146 non-smooth
```
