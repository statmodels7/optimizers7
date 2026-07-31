# S7 Class for Optimisation Algorithms

An optimiser is an object carrying an algorithm and every setting that
algorithm obeys: its stopping rule, its budgets, what it reports, and
what it keeps. It is created once and may be reused, inspected, stored
beside a result and passed around; nothing about a run is hidden in the
call that started it.

## Usage

``` r
optimizer(
  name = character(0),
  criterion = NULL,
  maxit = integer(0),
  max_eval = integer(0),
  verbose = logical(0),
  refresh = integer(0),
  keep_trace = logical(0)
)
```

## Arguments

- name:

  A short character name, used when printing and reporting.

- criterion:

  The stopping rule, a
  [`criterion`](https://statmodels7.github.io/optimizers7/reference/criterion.md)
  object.

- maxit:

  Maximum number of iterations.

- max_eval:

  Maximum number of objective evaluations. A budget on work rather than
  on progress: a line search can spend many evaluations in one
  iteration, and a runaway one is invisible to `maxit`.

- verbose:

  Logical; whether to report progress.

- refresh:

  Report every `refresh` iterations. `0` reports only the final summary.

- keep_trace:

  Logical; whether to store the iteration path in the result.

## Value

An S7 object of class `optimizer`.

## Details

The class is abstract. Each algorithm is a subclass with its own
constructor —
[`gd`](https://statmodels7.github.io/optimizers7/reference/gd.md), and
in time
[`newton()`](https://statmodels7.github.io/optimizers7/reference/Newton.md),
[`bfgs()`](https://statmodels7.github.io/optimizers7/reference/Bfgs.md)
and the rest — adding only the settings that are genuinely its own.

These properties are shared by every algorithm without exception, which
is what lets a caller swap one optimiser for another without changing
anything else.

## See also

[`minimize`](https://statmodels7.github.io/optimizers7/reference/minimize.md),
[`gd`](https://statmodels7.github.io/optimizers7/reference/gd.md)

## Examples

``` r
# Abstract: use one of the constructors.
try(optimizer(name = "mine"))
#> Error in new_object(S7_object(), name = name, criterion = criterion, maxit = maxit,  : 
#>   Can't construct an object from abstract class <optimizer>

# Every algorithm carries the same settings, which is what lets one be
# swapped for another without changing anything else.
names(S7::props(bfgs()))
#>  [1] "name"        "criterion"   "maxit"       "max_eval"    "verbose"    
#>  [6] "refresh"     "keep_trace"  "step"        "line_search" "curv_tol"   
#> [11] "max_skip"   
bfgs()@maxit
#> [1] 500
```
