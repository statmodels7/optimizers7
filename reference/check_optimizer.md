# Check What an Optimizer Reports

Puts an optimizer through twelve checks on what it *reports*, prints the
verdict on each, and then runs it over the standard test problems and
reports the gap that it reached on each one. It is intended for a caller
who adds a method. Every shipped method except
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
passes all twelve checks;
[`brent()`](https://statmodels7.github.io/optimizers7/reference/brent.md)
minimizes a function of one variable and cannot run the battery, whose
problems have more parameters.

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

  A list of problems in the shape that
  [`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
  returns. Defaults to every problem of
  [`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md).
  Pass a subset to keep the report short: `test_problems("sphere")`.

- verbose:

  `TRUE` (default) prints the twelve verdicts and the battery table;
  `FALSE` returns the same information silently.

- tol:

  Tolerance for the checks that compare numbers, a single positive
  number, default `1e-6`. It applies to check 1 and to the comparison of
  the two values in check 9. Check 2 uses a fixed `1e-3` on the
  gradient, check 9 compares the two points with a fixed relative
  tolerance of `1e-5`, check 12 uses a fixed `1e-3` on the distance to
  the solution, and check 7 tests the bounds strictly, without a
  tolerance.

## Value

Invisibly, a named list of two:

- `checks`:

  a named logical vector of length 12, one entry per numbered check,
  named by what the check asserts.

- `battery`:

  a data frame with one row per problem and the columns `problem`,
  `value`, `gap`, `converged`, `evaluations` and `note`. A problem that
  the optimizer could not run keeps its name in `problem`, has `NA` in
  the other numeric and logical columns, and has the error message in
  `note`.

## Requirements and performance

The twelve numbered checks are statements that an optimizer must satisfy
whatever its performance. The performance comes afterwards, as a table
of gaps without a verdict. The two are kept apart because a slow method
can satisfy every requirement: gradient descent does not solve
Rosenbrock in five hundred iterations, which is a documented property of
the method.

The only performance requirement among the numbered checks is check 12,
the minimization of a quadratic, which every correct method passes.

## The twelve checks

1.  `value` is the objective at `par`. A method that reports the value
    at a point that it has since left passes every test that uses only
    the value.

2.  the reported gradient is the gradient at `par`. The check is run
    only for an optimizer that declares `"gradient"` through
    [`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md).
    [`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
    reports an aggregate subgradient and does not declare a gradient, so
    the check does not apply to it.

3.  `converged` follows the stopping rule and is never inferred from the
    run having ended. The check starves the optimizer of iterations, and
    a run cut off after one iteration must not report success.

4.  budgets are respected: `iterations` never exceeds `maxit`.

5.  evaluations are counted: the reported number of evaluations is not
    zero.

6.  the trace, when kept, is a data frame whose key column (`iteration`,
    or `start` for a
    [`multistart()`](https://statmodels7.github.io/optimizers7/reference/multistart.md)
    trace) starts at one and increases.

7.  bounds are respected strictly. A probability of exactly 1 is not a
    probability inside \\(0, 1)\\, and the next computation of the
    caller often divides by it.

8.  the run repeats. A deterministic method gives the same answer twice;
    a stochastic one gives it again from the seed that it recorded,
    which tests the recording as well as the repeatability.

9.  [`maximize()`](https://statmodels7.github.io/optimizers7/reference/maximize.md)
    is
    [`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
    of the negative.

10. a stopping rule that the optimizer cannot evaluate is rejected when
    the run starts, and is not accepted and left never to fire.

11. a starting point where the objective is not finite raises an error,
    so that a run does not return `NaN` silently.

12. it minimizes a quadratic.

An optimizer that returns an invented point with `converged = TRUE`
fails six of the twelve checks: 1, 2, 3, 10, 11 and 12. The six that it
passes concern bookkeeping.

## The problem battery

The table reports the gap between the value reached and the known
minimum, as information. A large gap on `rastrigin` or `himmelblau`
means the method found a different local minimum, which for a local
method is correct behavior, and the `note` column labels those two. A
large gap on `abs_sum` means that the method stopped at a kink, the case
for which
[`bundle()`](https://statmodels7.github.io/optimizers7/reference/bundle.md)
and the derivative-free methods are intended.

## See also

[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
for the battery,
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
and
[`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
for the two functions that a new method has to call,
[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md)
for the declaration read by check 2.

## Examples

``` r
check_optimizer(bfgs())
#> Checking optimizer: BFGS
#>   [ 1] value agrees with par:                       [PASSED]
#>   [ 2] gradient agrees with par:                    [PASSED]
#>   [ 3] convergence is not assumed:                  [PASSED]
#>   [ 4] budgets are respected:                       [PASSED]
#>   [ 5] evaluations are counted:                     [PASSED]
#>   [ 6] trace is well formed:                        [PASSED]
#>   [ 7] bounds are respected strictly:               [PASSED]
#>   [ 8] the run repeats:                             [PASSED]
#>   [ 9] maximize mirrors minimize:                   [PASSED]
#>   [10] a rule that cannot be evaluated is rejected: [PASSED]
#>   [11] a bad starting point is an error:            [PASSED]
#>   [12] it minimizes a quadratic:                    [PASSED]
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

# For a method that computes no gradient, the gradient check does not
# apply, and the other checks are the same.
res <- check_optimizer(nelder_mead(), problems = test_problems("sphere"))
#> Checking optimizer: nelder-mead
#>   [ 1] value agrees with par:                       [PASSED]
#>   [ 2] gradient agrees with par:                    [PASSED]
#>   [ 3] convergence is not assumed:                  [PASSED]
#>   [ 4] budgets are respected:                       [PASSED]
#>   [ 5] evaluations are counted:                     [PASSED]
#>   [ 6] trace is well formed:                        [PASSED]
#>   [ 7] bounds are respected strictly:               [PASSED]
#>   [ 8] the run repeats:                             [PASSED]
#>   [ 9] maximize mirrors minimize:                   [PASSED]
#>   [10] a rule that cannot be evaluated is rejected: [PASSED]
#>   [11] a bad starting point is an error:            [PASSED]
#>   [12] it minimizes a quadratic:                    [PASSED]
#> 
#>   All checks passed.
#> 
#>   battery (gap from the known minimum; information, not a verdict)
#>     sphere       gap  2.14e-17  conv      276 evals  
all(res$checks)
#> [1] TRUE

# The battery is returned as a table.
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
