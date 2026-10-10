# Extending optimizers7

This vignette describes how to add a stopping rule, a kind of objective
and an algorithm to the package, and which checks the package applies to
the added pieces.

The three components are separate, so that replacing one of them
requires no change to the other two.

## A user-defined stopping rule

The alternative to an object would be an argument that takes a string,
with a `switch` inside every algorithm. That arrangement fixes the set
of meanings of *finished* when the package is written, and the meaning
of convergence determines what the reported answer is.

A criterion is one class and one method. The method receives the state
of the iteration and returns a single logical:

``` r

CritTiny <- S7::new_class("CritTiny", parent = criterion,
                          properties = list(tol = S7::class_numeric))

S7::method(crit_met, CritTiny) <- function(criterion, state) {
  state$f_new < criterion@tol
}

tiny <- CritTiny(label = "objective below 1e-8", tol = 1e-8)
```

Every algorithm in the package then accepts it:

``` r

f  <- function(p) sum((p - c(1, 2))^2)
gr <- function(p) 2 * (p - c(1, 2))

minimize(bfgs(criterion = tiny), f, c(0, 0), gr = gr)@criterion_met
#> [1] "objective below 1e-8"
```

`state` carries `iter`, `f_new`, `f_old`, `x_new`, `x_old`, `gradient`
and `stationarity`. A method supplies `gradient` or `stationarity`
according to its kind and `NULL` for the other, and `f_old` and `x_old`
are `NULL` when there is no previous value; the help page of
[`crit_met()`](https://statmodels7.github.io/optimizers7/reference/crit_met.md)
gives the details.

### Declaring what a rule needs

A rule that reads the gradient, used with a method that computes none,
would test `NULL` at every iteration and never fire; the run would end
on its iteration budget and report a failure unrelated to the state of
the run. A rule therefore declares what it reads, and an optimizer that
cannot supply it rejects the rule with an error that names the missing
quantity:

``` r

minimize(nelder_mead(criterion = crit_grad()), f, c(0, 0))
#> Error:
#> ! The stopping rule needs gradient, which nelder-mead does not provide.
#>   Choose a criterion this optimizer can evaluate, or a method that provides it.
```

A rule declares the components that it needs through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md);
a rule that reads only the objective, as `CritTiny` does, declares
nothing:

``` r

S7::method(crit_needs, CritTiny) <- function(criterion) character()
crit_needs(crit_grad())
#> [1] "gradient"
crit_needs(crit_stationary())
#> [1] "stationarity"
```

Rules combine, and the combination is itself a rule, so combinations
nest:

``` r

crit_any(tiny, crit_grad(1e-10))
#> <criterion> objective below 1e-8 or gradient (max-norm) < 1e-10
```

## The objective

[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
dispatches on the optimizer, so an algorithm is written once. The
objective is normalized separately, by
[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md),
which dispatches on the objective. Dispatching
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
on both would need one method for each pair of an algorithm and a kind
of objective.

### An ordinary function

This is the common case. Every evaluation is a callback into R.

``` r

minimize(bfgs(), f, c(0, 0), gr = gr)@counts
#> f g h 
#> 3 3 0
```

With no gradient supplied, one is obtained by central differences, and
the result records it:

``` r

minimize(bfgs(), f, c(0, 0))@message
#> [1] "gradient obtained by finite differences"
```

### An objective that resamples itself

Subsampling is a property of the objective, not of an algorithm, and it
needs no support from the package. An objective that draws a new
minibatch on every call is an ordinary R function. Here `batch` and
`batch_gr` draw different minibatches, so the gradient check of
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
gives a warning, which `options(optimizers7.check_gradient = FALSE)`
switches off:

``` r

y <- rnorm(2000, mean = 3)
m <- 100

batch    <- function(par) { i <- sample.int(2000, m); sum((y[i] - par)^2) / 2 }
batch_gr <- function(par) { i <- sample.int(2000, m); -sum(y[i] - par) }

minimize(adam(alpha = 0.05, decay = 0.005, maxit = 3000),
         batch, par = 0, gr = batch_gr)@par
#> Warning: 'gr' does not appear to be the gradient of 'fn': along the gradient direction at 'par',
#>   'fn' changes at rate -2.4e+08 where 'gr' predicts 291. Check that the two
#>   compute the same model. options(optimizers7.check_gradient = FALSE) turns this check off.
#> [1] 2.987537
mean(y)
#> [1] 2.986045
```

The package has no objective class that declares itself a sum over
observations for
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
to subset. Such a class would require a second kind of objective, a rule
for the stopping rules that it permits, a way of reporting which one is
in force, and a branch in the compiled loop, for a task that the caller
can do in one line, since the caller defines what an observation is. An
optimizer with its own notion of an observation would not be a general
optimizer, so the batch is drawn in the layer that defines observations.

The resampling belongs *inside* the objective, not around the run.
Calling `minimize(adam(maxit = 1), ...)` in a loop resets $`m`$ and
$`v`$ to zero and the bias correction to $`t = 1`$ at every call. Each
call then takes a first step of length $`\alpha`$ in each coordinate,
and the accumulated moments, on which the method depends, are discarded
at each call.

### A user-defined kind of objective

The package ships one kind of objective, and
[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md)
is a generic so that there can be others. A caller with a fitted-model
object, for example, registers one method and every algorithm in the
package accepts the object:

``` r

S7::method(as_objective, MyModel) <- function(fn, gr = NULL, he = NULL, ...) {
  list(kind = "r",
       fn = function(par) -loglik(fn, par),
       gr = function(par) -score(fn, par),
       he = NULL,
       has_gradient = TRUE, has_hessian = FALSE)
}
```

The package does not provide a compiled objective, a pair of external
pointers through which the iteration would never return to R, for two
reasons.

First, such an objective cannot carry data. The pointer type would be
`double(*)(const arma::vec&)`, a bare function pointer with no closure
and no user-data argument, so a statistical objective would have to keep
its `y` and its `X` in C++ globals. That is not reentrant, it allows one
model at a time, and it rules out running several starts in one process.

Second, it gives no general speed gain. The gradient and the Hessian of
a model are assembled from `crossprod(Z, s)` and `crossprod(Z, W * Z)`,
and R and Armadillo call the same BLAS for these, so the dominant
operation runs the same code in both languages. Compiled kernels pay
where the arithmetic is elementwise and irregular, which is where
distributions7 uses them, in the derivative expressions of each family
up to the fifth order.

## A user-defined algorithm

An optimizer is a subclass of `optimizer` that carries its own settings,
and a method of
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
The example is a method that the package does not ship, the heavy-ball
method (gradient descent with momentum), which adds a fraction of the
previous *displacement* to each step:

``` math
x_{k+1} = x_k - \alpha g_k + \beta (x_k - x_{k-1}).
```

The momentum term carries the iterate along the floor of a valley
instead of across it, and it is the precursor of the moment estimate in
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md).

``` r

HeavyBall <- S7::new_class("HeavyBall", parent = optimizer,
                           properties = list(alpha = S7::class_numeric,
                                             beta  = S7::class_numeric))

heavy_ball <- function(criterion = crit_grad(1e-8), alpha = 1e-3, beta = 0.9,
                       maxit = 2000, max_eval = 20000, verbose = FALSE,
                       refresh = 50, keep_trace = FALSE) {
  HeavyBall(name = "heavy ball", criterion = criterion,
            maxit = maxit, max_eval = max_eval, verbose = verbose,
            refresh = refresh, keep_trace = keep_trace,
            alpha = alpha, beta = beta)
}

hb_run <- function(o, fn, par, gr) {
  x <- as.numeric(par); x_prev <- x
  nf <- 0L; ng <- 0L; converged <- FALSE; it <- 0L

  for (it in seq_len(o@maxit)) {
    g <- gr(x); ng <- ng + 1L
    x_new <- x - o@alpha * g + o@beta * (x - x_prev)
    # The safeguard, and it is not optional: momentum can carry the iterate
    # somewhere the objective is undefined, and one NaN contaminates every
    # iterate after it.
    if (!all(is.finite(x_new))) break
    x_prev <- x; x <- x_new

    state <- list(iter = it, f_new = fn(x), f_old = fn(x_prev), x_new = x,
                  x_old = x_prev, gradient = gr(x))
    nf <- nf + 2L; ng <- ng + 1L
    if (crit_met(o@criterion, state)) { converged <- TRUE; break }
  }

  optimizer_result(
    par = x, value = fn(x), gradient = gr(x),
    counts = c(f = nf + 1L, g = ng + 1L, h = 0L), iterations = it,
    converged = converged,
    criterion_met = if (converged) o@criterion@label
                    else "iteration budget reached",
    message = "", trace = NULL, optimizer = o, elapsed = NA_real_
  )
}

S7::method(minimize, HeavyBall) <-
  function(optimizer, fn, par, gr = NULL, he = NULL, lower = -Inf, upper = Inf, ...) {
    if (is.null(gr)) stop("heavy_ball() needs a gradient.", call. = FALSE)
    hb_run(optimizer, fn, par, gr)
  }
```

The method sets `converged` itself: it is `TRUE` when the stopping rule
fired and `FALSE` otherwise, including when the run ended on its
iteration budget.

On a quadratic:

``` r

r <- minimize(heavy_ball(alpha = 0.05), f, c(0, 0), gr = gr)
c(par = r@par, iterations = r@iterations, converged = r@converged)
#>       par1       par2 iterations  converged 
#>          1          2        327          1
```

## Checking it

[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
separates two questions: whether the optimizer satisfies the twelve
numbered checks, and how **powerful** it is. A slow method is not a
defective one: gradient descent does not solve Rosenbrock in five
hundred iterations, and that is a property of the method.

The numbered checks are therefore statements that every optimizer must
satisfy, whatever its performance, and the battery below them is a table
of gaps, reported as information and not as a verdict:

``` r

res <- check_optimizer(heavy_ball(alpha = 0.01))
#> Checking optimizer: heavy ball
#>   [ 1] value agrees with par:            [PASSED]
#>   [ 2] gradient agrees with par:         [PASSED]
#>   [ 3] convergence is not assumed:       [PASSED]
#>   [ 4] budgets are respected:            [PASSED]
#>   [ 5] evaluations are counted:          [PASSED]
#>   [ 6] trace is well formed:             [PASSED]
#>   [ 7] bounds are respected strictly:    [FAILED]
#>   [ 8] the run repeats:                  [PASSED]
#>   [ 9] maximize mirrors minimize:        [PASSED]
#>   [10] an unevaluable rule is rejected:  [FAILED]
#>   [11] a bad starting point is an error: [PASSED]
#>   [12] it minimizes a quadratic:         [PASSED]
#> 
#>   2 check(s) FAILED: bounds are respected strictly; an unevaluable rule is rejected
#> 
#>   battery (gap from the known minimum; information, not a verdict)
#>     sphere       gap  8.47e-18  conv      677 evals  
#>     rosenbrock   gap       Inf  -          15 evals  
#>     booth        gap  2.07e-17  conv      743 evals  
#>     beale        gap  3.03e-17  conv      755 evals  
#>     powell       gap       Inf  -          13 evals  
#>     himmelblau   gap  8.60e-19  conv      861 evals  multimodal
#>     rastrigin    gap  4.46e+00  -        4001 evals  multimodal
#>     abs_sum      gap  3.39e-02  -        4001 evals  non-smooth
```

The battery shows a limitation of the method: with the fixed step length
of 0.01 the iterate diverges on `rosenbrock` and `powell`, whose
curvature varies by orders of magnitude, and the run does not converge
on `rastrigin` and `abs_sum`. This is a property of heavy ball with a
fixed step.
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md),
which the package ships, has the same form with the step length
estimated from the last secant pair, and it converges on all eight
problems.

The two failed checks concern the method and not the algorithm: the
method receives two arguments and does nothing with them, and it does
not validate its stopping rule.

**It ignores `lower` and `upper`.** The signature includes them, since
S7 warns when the formals of a method omit a named argument of the
generic, and an argument that the author did not plan for is easily left
unused. A caller then passes a box, gets no error, and receives an
answer outside it.

Removing the box is a few lines, because
[`check_bounds()`](https://statmodels7.github.io/optimizers7/reference/check_bounds.md)
normalizes the two vectors into one pair per coordinate and
[`bounded_transform()`](https://statmodels7.github.io/optimizers7/reference/bounded_transform.md)
does the rest:

``` r

suppressMessages(
S7::method(minimize, HeavyBall) <-
  function(optimizer, fn, par, gr = NULL, he = NULL, lower = -Inf, upper = Inf, ...) {
    if (is.null(gr)) stop("heavy_ball() needs a gradient.", call. = FALSE)
    check_criterion(optimizer)                                     # (2)

    bx <- check_bounds(lower, upper, par)                       # (1)
    if (!length(bx)) return(hb_run(optimizer, fn, par, gr))

    # Optimize in eta, where there is no box, and map back to report. The
    # Jacobian is diagonal, so the chain rule is one multiplication.
    hh <- function(e) vapply(seq_along(e),
            function(j) bounded_transform(bx[[j]], e[j])$h, 0)
    dd <- function(e) vapply(seq_along(e),
            function(j) bounded_transform(bx[[j]], e[j])$d1, 0)
    eta0 <- vapply(seq_along(par),
                   function(j) bounded_forward(bx[[j]], par[j]), 0)

    res <- hb_run(optimizer,
                  function(e) fn(hh(e)),
                  eta0,
                  function(e) gr(hh(e)) * dd(e))

    res@gradient <- res@gradient / dd(res@par)   # both on the user's scale
    res@par <- hh(res@par)
    res
  }
)
```

**It does not validate the stopping rule.**
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
(marked (2)) is the line that rejects a rule that the method cannot
evaluate, instead of accepting one that would test `NULL` at every
iteration. The quantities that a method supplies are declared by
[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md),
whose default method returns `"gradient"`. This method computes the
gradient, so the default is correct for it. A derivative-free method
would need one line:

``` r

S7::method(optimizer_provides, HeavyBall) <- function(optimizer)
  "stationarity"
```

With both repairs, every check passes:

``` r

all(check_optimizer(heavy_ball(alpha = 0.01),
                    problems = test_problems("sphere"),
                    verbose = FALSE)$checks)
#> [1] TRUE
```

The two repairs add the handling of bounds and the validation of the
stopping rule, which the package requires of every optimizer. Neither
omission affects a run without bounds and with a rule that the method
can evaluate, and the numbered checks are written to detect both.

## The test problems

[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
is exported because a battery is useful for writing any optimizer, not
only those of this package:

``` r

p <- test_problems("beale")[[1]]
str(p, max.level = 1)
#> List of 9
#>  $ name      : chr "beale"
#>  $ fn        :function (p)  
#>  $ gr        :function (p)  
#>  $ he        :function (p)  
#>  $ par       : num [1:2] 1 1
#>  $ solution  : num [1:2] 3 0.5
#>  $ value     : num 0
#>  $ multimodal: logi FALSE
#>  $ smooth    : logi TRUE

minimize(bfgs(), p$fn, p$par, gr = p$gr)@par
#> [1] 3.0 0.5
p$solution
#> [1] 3.0 0.5
```

Two of them are marked `multimodal`: on them a local method that reaches
a different minimum is behaving correctly. One is marked non-`smooth`:
on it a method that assumes a derivative reaches the solution and cannot
certify it:

``` r

sapply(test_problems(), function(p) c(multimodal = p$multimodal,
                                      smooth = p$smooth))
#>            sphere rosenbrock booth beale powell himmelblau rastrigin abs_sum
#> multimodal  FALSE      FALSE FALSE FALSE  FALSE       TRUE      TRUE   FALSE
#> smooth       TRUE       TRUE  TRUE  TRUE   TRUE       TRUE      TRUE   FALSE
```
