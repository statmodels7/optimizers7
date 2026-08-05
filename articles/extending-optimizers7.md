# Extending optimizers7

A toolkit that only worked for the algorithms it ships with would have
solved nothing. This vignette covers the other half: what can be added,
and what the package guarantees about the added pieces.

Three components can be replaced: the stopping rule, the objective, and
the algorithm. They are separate on purpose, so that replacing one costs
nothing in the other two.

## A user-defined stopping rule

The alternative to an object here would be an argument taking a string
and a `switch` inside every algorithm. That arrangement fixes, at the
moment the package is written, what everyone downstream is allowed to
mean by *finished*, and convergence is the definition of the answer
rather than a detail.

A criterion is one class and one method. The method receives the whole
state of the iteration and returns a single logical:

``` r

CritTiny <- S7::new_class("CritTiny", parent = criterion,
                          properties = list(tol = S7::class_numeric))

S7::method(crit_met, CritTiny) <- function(criterion, state) {
  state$f_new < criterion@tol
}

tiny <- CritTiny(label = "objective below 1e-8", tol = 1e-8)
```

It is then usable by every algorithm in the package, with no further
ceremony:

``` r

f  <- function(p) sum((p - c(1, 2))^2)
gr <- function(p) 2 * (p - c(1, 2))

minimize(bfgs(criterion = tiny), f, c(0, 0), gr = gr)@criterion_met
#> [1] "objective below 1e-8"
```

`state` carries `iter`, `f_new`, `f_old`, `x_new`, `x_old`, `gradient`
and `stationarity`. The last two are not always there, and that matters.

### Declaring what a rule needs

A rule that reads the gradient, handed to a method that computes none,
would sit testing `NULL` at every iteration and quietly never fire; the
run would end on its iteration budget and report a reason nowhere near
the truth. So a rule declares what it reads, and an optimizer that
cannot supply it refuses the rule by name:

``` r

minimize(nelder_mead(criterion = crit_grad()), f, c(0, 0))
#> Error:
#> ! The stopping rule needs gradient, which nelder-mead does not provide.
#>   Choose a criterion this optimizer can evaluate, or a method that provides it.
```

A rule that reads something optional declares it through
[`crit_needs()`](https://statmodels7.github.io/optimizers7/reference/crit_needs.md):

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
which dispatches on the objective. Each generic dispatches where there
is real variation; dispatching
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md)
on both would need one method per algorithm per shape.

### An ordinary function

The common case, and the one that costs the most: every evaluation is a
callback into R.

``` r

minimize(bfgs(), f, c(0, 0), gr = gr)@counts
#> f g h 
#> 3 2 0
```

With no gradient supplied, one is obtained by central differences, and
the result says so, so a run is never silently less exact than it
appears:

``` r

minimize(bfgs(), f, c(0, 0))@message
#> [1] "gradient obtained by finite differences"
```

### An objective that resamples itself

Subsampling is a property of the objective, not of an algorithm, and the
sharpest way to say so is that it needs no support at all. An objective
that draws a fresh minibatch on every call is an ordinary R function:

``` r

y <- rnorm(2000, mean = 3)
m <- 100

batch    <- function(par) { i <- sample.int(2000, m); sum((y[i] - par)^2) / 2 }
batch_gr <- function(par) { i <- sample.int(2000, m); -sum(y[i] - par) }

minimize(adam(alpha = 0.05, decay = 0.005, maxit = 3000),
         batch, par = 0, gr = batch_gr)@par
#> Warning: 'gr' does not appear to be the gradient of 'fn': along the gradient direction at 'par',
#>   'fn' changes at rate -3316475 where 'gr' predicts 291. Check that the two
#>   compute the same model. options(optimizers7.check_gradient = FALSE) turns this check off.
#> [1] 2.986939
mean(y)
#> [1] 2.986045
```

An earlier version of the package had a `finite_sum()` class, by which
an objective declared that it was a sum over observations so that
[`adam()`](https://statmodels7.github.io/optimizers7/reference/adam.md)
could ask it for a subset. It was removed. Supporting it meant a second
kind of objective, a rule for which stopping rules that objective
permitted, a way of reporting which was in force, and a branch in the
compiled loop. All of that served something the caller can do in one
line and do better, since the caller is the one who knows what an
observation is. An optimizer that knows about observations has stopped
being a general one; the batch is drawn where the knowledge is.

One detail matters when writing such an objective. The resampling
belongs *inside* the objective, not around the run. Calling
`minimize(adam(maxit = 1), ...)` in a loop resets and to zero and the
bias correction to at every call. Each call then takes a first step of
length , and the accumulated moments, which are the whole of the method,
are thrown away each time.

### A user-defined kind of objective

There is one shipped shape, and
[`as_objective()`](https://statmodels7.github.io/optimizers7/reference/as_objective.md)
is a generic so that there can be others. A caller holding a
fitted-model object, say, registers one method and every algorithm in
the package accepts it:

``` r

S7::method(as_objective, MyModel) <- function(fn, gr = NULL, he = NULL, ...) {
  list(kind = "r",
       fn = function(par) -loglik(fn, par),
       gr = function(par) -score(fn, par),
       he = NULL,
       has_gradient = TRUE, has_hessian = FALSE)
}
```

An earlier version also shipped a second shape, `cpp_objective()`,
taking a pair of external pointers so that the iteration never returned
to R. It was removed for two reasons.

First, it could not carry data. The pointer type was
`double(*)(const arma::vec&)`, a bare function pointer with no closure
and no user-data argument, so any real statistical objective had to keep
its `y` and its `X` in C++ globals. That is not reentrant, it means one
model at a time, and it rules out running several starts in one process.

Second, it was not faster where it matters. Measured on a Gamma
regression against the same objective written in vectorized R, the
compiled version was *slower* below about twenty thousand observations
and never better than about twice as fast above it. The reason is
structural: assembling a gradient or a Hessian from model terms is
`crossprod(Z, s)` and `crossprod(Z, W * Z)`, and R’s matrix arithmetic
and Armadillo’s call the same BLAS, so on the operation that dominates a
real fit the two languages run the same code. Compiled kernels pay where
the arithmetic is elementwise and irregular, which is where puts them,
in the fourth-order derivative expressions of each family.

## A user-defined algorithm

An optimizer is a subclass of `optimizer` carrying its own settings,
plus a method on
[`minimize()`](https://statmodels7.github.io/optimizers7/reference/minimize.md).
Here is one the package does not ship — the heavy-ball method, gradient
descent with momentum, which adds a fraction of the previous
*displacement* to each step:

``` math
x_{k+1} = x_k - \alpha g_k + \beta (x_k - x_{k-1}).
```

The momentum term is what carries the iterate along a valley floor
instead of across it, and it is the ancestor of the moment estimate in
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

Note the thing the method is responsible for and that nothing can do for
it: `converged` is set from whether the rule fired and **never** from
the run having ended.

It works:

``` r

r <- minimize(heavy_ball(alpha = 0.05), f, c(0, 0), gr = gr)
c(par = r@par, iterations = r@iterations, converged = r@converged)
#>       par1       par2 iterations  converged 
#>          1          2        327          1
```

## Checking it

[`check_optimizer()`](https://statmodels7.github.io/optimizers7/reference/check_optimizer.md)
separates two questions that are easy to confuse: whether the optimizer
keeps its **contract**, and how **powerful** it is. Conflating them
would make the function useless, since gradient descent does not solve
Rosenbrock in five hundred iterations and is not broken but slow.

The numbered checks are therefore statements every optimizer must
satisfy however weak it is, and the battery underneath is a table of
gaps, reported as information rather than as a verdict:

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
#>   [10] an unevaluable rule is refused:   [FAILED]
#>   [11] a bad starting point is an error: [PASSED]
#>   [12] it minimizes a quadratic:         [PASSED]
#> 
#>   2 check(s) FAILED: bounds are respected strictly; an unevaluable rule is refused
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

The battery shows a real limitation: a fixed step length makes this
method crawl on problems whose curvature varies by orders of magnitude,
and that is a fair description of heavy ball rather than a defect in it.
[`bb()`](https://statmodels7.github.io/optimizers7/reference/bb.md),
which the package does ship, is the same shape of method with the step
length estimated from the last secant pair instead of fixed, and it
solves those same problems to machine precision.

The two failed checks are the more useful half of the output. Nothing
above is wrong with the *algorithm*; the method accepted two arguments
and did nothing with them.

**It ignored `lower` and `upper`.** The signature has to include them,
since S7 requires a method’s formals to contain every named argument of
the generic, and an argument the author did not plan for is easily left
sitting there. A caller then passes a box, gets no error, and receives
an answer outside it.

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

**It never validated the stopping rule.**
[`check_criterion()`](https://statmodels7.github.io/optimizers7/reference/check_criterion.md)
— marked (2) — is the one line that refuses a rule the method cannot
evaluate, instead of accepting one that would test `NULL` for ever. What
a method *can* supply is declared by
[`optimizer_provides()`](https://statmodels7.github.io/optimizers7/reference/optimizer_provides.md),
whose default claims a gradient and an objective. This one has both, so
the default is honest; a derivative-free method would need one line:

``` r

S7::method(optimizer_provides, HeavyBall) <- function(optimizer)
  c("objective", "stationarity")
```

With both promises kept, the contract holds:

``` r

all(check_optimizer(heavy_ball(alpha = 0.01),
                    problems = test_problems("sphere"),
                    verbose = FALSE)$checks)
#> [1] TRUE
```

The algorithm was right from the first line. What was missing was two
promises the package makes on every optimizer’s behalf, and neither
omission would have shown itself in any run that happened to succeed;
the numbered checks exist to expose exactly this kind of gap.

## The test problems

[`test_problems()`](https://statmodels7.github.io/optimizers7/reference/test_problems.md)
is exported because a battery is useful to anyone writing an optimizer,
not only to this package:

``` r

p <- test_problems("beale")[[1]]
str(p, max.level = 1)
#> List of 8
#>  $ name      : chr "beale"
#>  $ fn        :function (p)  
#>  $ gr        :function (p)  
#>  $ par       : num [1:2] 1 1
#>  $ solution  : num [1:2] 3 0.5
#>  $ value     : num 0
#>  $ multimodal: logi FALSE
#>  $ smooth    : logi TRUE

minimize(bfgs(), p$fn, p$par, gr = p$gr)@par
#> [1] 3.0000005 0.5000001
p$solution
#> [1] 3.0 0.5
```

Two of them are marked `multimodal`, where a local method reaching a
different minimum is behaving correctly rather than failing. One is
marked non-`smooth`, where a method that assumes a derivative exists
will walk to the answer and then be unable to certify it:

``` r

sapply(test_problems(), function(p) c(multimodal = p$multimodal,
                                      smooth = p$smooth))
#>            sphere rosenbrock booth beale powell himmelblau rastrigin abs_sum
#> multimodal  FALSE      FALSE FALSE FALSE  FALSE       TRUE      TRUE   FALSE
#> smooth       TRUE       TRUE  TRUE  TRUE   TRUE       TRUE      TRUE   FALSE
```
