#' @include generics.R
#' @include result_class.R
NULL

# The two pieces every algorithm's method needs, written once here so that a new
# algorithm is a C++ loop and a constructor and nothing else.


#' Prepare an Objective for an Algorithm
#'
#' @description
#' Normalizes the objective, checks that this optimizer can evaluate the
#' stopping rule, and checks the starting value.
#'
#' @details
#' A rule that needs a gradient, used with a method that computes none,
#' would test `NULL` at every iteration and never fire, so the run would end
#' on its iteration budget and report a failure unrelated to the state of the
#' run. The check rejects such a rule when the run starts, with an error that
#' names the missing quantity and the method.
#'
#' @param optimizer The [optimizer()].
#' @param fn,gr,he The objective and its optional derivatives, as supplied.
#' @param par The starting value.
#'
#' @return The objective handle from [as_objective()].
#'
#' @keywords internal
prepare_objective <- function(optimizer, fn, par, gr = NULL, he = NULL) {
  if (!is.numeric(par) || !length(par) || anyNA(par)) {
    stop("'par' must be a numeric vector of starting values.", call. = FALSE)
  }
  spec <- as_objective(fn, gr, he)

  check_criterion(optimizer)
  spec
}


#' @title Reject a Stopping Rule That an Optimizer Cannot Evaluate
#'
#' @description
#' Compares the quantities that the criterion of the optimizer reads with the
#' quantities that the optimizer supplies, and signals an error if the rule
#' needs a quantity that the optimizer does not supply.
#'
#' @param optimizer An [optimizer()].
#'
#' @details
#' Call this at the top of a user-defined [minimize()] method. A rule that
#' needs a gradient, used with a method that computes none, would test `NULL`
#' at every iteration and never fire; the run would then end on its iteration
#' budget and report a failure unrelated to the state of the run.
#' [check_optimizer()] tests that a method rejects such a rule.
#'
#' The quantities that an optimizer supplies are declared by
#' [optimizer_provides()], whose default method returns `"gradient"`.
#' Override it when the method computes no gradient.
#'
#' @return Invisibly `TRUE`; raises an error otherwise.
#'
#' @examples
#' check_criterion(bfgs())
#' try(check_criterion(nelder_mead(criterion = crit_grad())))
#'
#' @seealso [optimizer_provides()], [crit_needs()]
#' @export
check_criterion <- function(optimizer) {
  needs <- crit_needs(optimizer@criterion)
  provides <- optimizer_provides(optimizer)
  missing <- setdiff(needs, provides)
  if (length(missing)) {
    msg <- sprintf(
      paste0("The stopping rule needs %s, which %s does not provide.\n",
             "  Choose a criterion that this optimizer can evaluate, or a method ",
             "that provides it."),
      paste(missing, collapse = ", "), optimizer@name)
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}


#' @title What an Optimizer Can Offer a Stopping Rule
#'
#' @description
#' Returns the names of the `state` components that an optimizer supplies,
#' so that [check_criterion()] can reject a rule that the optimizer cannot
#' evaluate.
#'
#' @details
#' The gradient-based methods provide `"gradient"`. [nelder_mead()],
#' [compass()], [brent()], [sa()] and [bundle()] provide `"stationarity"`
#' instead: [nelder_mead()] offers the largest distance from the best vertex
#' to another vertex of the simplex, [compass()] the poll size, [brent()] the
#' half-width of its bracket, [sa()] Corana's termination measure and
#' [bundle()] its optimality estimate. [bundle()] evaluates subgradients but
#' does not provide `"gradient"`, since a subgradient does not vanish at a
#' kink.
#'
#' Every optimizer evaluates the objective, so there is no name for it, and a
#' rule that reads the objective is accepted by every method. [prox_grad()]
#' returns `"gradient"`, and its gradient is the proximal gradient mapping,
#' which vanishes at a stationary point of the whole objective.
#'
#' The default method returns `"gradient"`, so an optimizer written
#' outside the package declares a gradient unless it has a method of its
#' own. The rejection works only if the declaration is accurate: a method
#' that declares a gradient that it does not compute accepts a rule that can
#' never fire.
#'
#' @param optimizer An [optimizer()].
#'
#' @return A character vector of `state` component names. The two that the
#'   shipped criteria read are `"gradient"` and `"stationarity"`.
#'
#' @examples
#' # The shipped methods give both answers.
#' vapply(list(bfgs(), newton(), cg(), bb(), gd(), adam()),
#'        optimizer_provides, "")
#' vapply(list(nelder_mead(), compass(), sa(), bundle()),
#'        optimizer_provides, "")
#'
#' # A wrapper returns the answer of the run that reports the result.
#' optimizer_provides(chain(sa(), bfgs()))
#' optimizer_provides(multistart(nelder_mead()))
#'
#' @seealso [check_criterion()] for the rejection that this generic feeds,
#'   [crit_needs()] for the other half of the comparison,
#'   [optimizer_bounded()] for the package's other declaration generic.
#' @export
#' @aliases optimizer_provides.optimizer
optimizer_provides <- S7::new_generic("optimizer_provides", "optimizer",
                                      function(optimizer) S7::S7_dispatch())

S7::method(optimizer_provides, optimizer) <- function(optimizer) "gradient"


#' @title Whether an Optimizer Takes Box Bounds
#'
#' @description
#' Returns `TRUE` when the optimizer applies the `lower` and `upper`
#' arguments of [minimize()], so that [check_optimizer()] tests them, and
#' `FALSE` for a method that takes its constraint in another way.
#'
#' @details
#' Every method in the package except [prox_grad()] removes bounds by
#' reparametrization and returns `TRUE`, which is the default. [prox_grad()]
#' receives a constraint inside the proximal operator, where it combines
#' with the non-smooth term, so bounds beside the objective would be a
#' second, conflicting way of imposing the same constraint.
#'
#' A wrapper returns the answer of what it contains. [chain()] is bounded
#' only when every stage is, since the bounds are passed to all of them,
#' and [multistart()] returns the answer of its inner optimizer.
#'
#' [check_optimizer()] reads the declaration before its seventh check, so a
#' method that returns `FALSE` is not failed for rejecting bounds.
#'
#' @param optimizer An [optimizer()].
#'
#' @return A single logical.
#'
#' @examples
#' optimizer_bounded(bfgs())
#'
#' pg <- prox_grad(prox = function(v, t) v, g = function(b) 0)
#' optimizer_bounded(pg)
#'
#' # A chain is bounded only if every stage is.
#' c(optimizer_bounded(chain(bfgs(), newton())),
#'   optimizer_bounded(chain(bfgs(), pg)))
#'
#' @seealso [optimizer_provides()] for the other declaration generic,
#'   [check_bounds()] for what a method that returns `TRUE` receives,
#'   [prox_grad()] for the one that returns `FALSE`.
#' @export
#' @aliases optimizer_bounded.optimizer
optimizer_bounded <- S7::new_generic("optimizer_bounded", "optimizer",
                                     function(optimizer) S7::S7_dispatch())

S7::method(optimizer_bounded, optimizer) <- function(optimizer) TRUE


#' Translate an Evaluation Budget for the Compiled Loop
#'
#' @description
#' The compiled loops hold the budget as an integer, and the default budget is
#' `Inf`, which sets no cap at all. This maps `Inf` to the largest
#' representable integer, which a run never reaches, and any finite value to
#' itself.
#'
#' @param x The `max_eval` property, a single positive number or `Inf`.
#' @return A single integer.
#' @keywords internal
budget_int <- function(x) {
  if (is.finite(x)) as.integer(x) else .Machine$integer.max
}


#' Warn When a Supplied Gradient Does Not Belong to the Objective
#'
#' @description
#' Compares the directional derivative of `fn` along the gradient direction
#' at `par`, obtained from one central difference, with the rate predicted
#' by the supplied `gr`. A gradient computed from a different model than the
#' objective, such as one that holds a predictor fixed where the objective
#' varies it or that reads an outdated copy of the data, makes the two
#' disagree by a large factor. Without this check the effect is a
#' line-search failure at the first iteration, with no indication of its
#' cause.
#'
#' @details
#' The check costs one call to `gr` and two to `fn`, runs once per call of
#' [minimize()], and is skipped when the comparison cannot be made: for an
#' objective that is not a function, a gradient that is not finite at `par`,
#' an objective that is not finite at the probe points, and a gradient too
#' small to be distinguished from the truncation error of the difference.
#' The last case arises when a caller starts at the optimum, where a warning
#' would be wrong.
#'
#' The tolerance is loose, a relative disagreement above one half, so that
#' the error of the difference or a subgradient of a non-smooth objective
#' does not trigger it. It is meant to detect a wrong function, not an error
#' in the eighth digit. A disagreement at the first step is checked again at
#' a step ten times smaller and, if it persists, at a step a hundred times
#' smaller, and the warning is given only if it persists at both: the step
#' is scaled by the length of `par`, so a large coordinate can make it too
#' long for a strongly curved direction, and the truncation error of the
#' difference shrinks with the step, whereas the disagreement of a wrong
#' gradient does not. Each repetition costs two more calls to `fn`.
#'
#' Setting `options(optimizers7.check_gradient = FALSE)` disables it.
#' [multistart()] disables it for all its inner runs, since the generic has
#' already checked the gradient at the starting value of the caller.
#'
#' @param fn The objective.
#' @param gr The supplied gradient.
#' @param par The starting point, numeric.
#'
#' @return Invisibly `TRUE`; warns on a gross mismatch.
#'
#' @keywords internal
check_gradient_consistency <- function(fn, gr, par) {
  if (!isTRUE(getOption("optimizers7.check_gradient", TRUE))) {
    return(invisible(TRUE))
  }
  if (!is.function(fn) || !is.function(gr) || !is.numeric(par) || !length(par)) {
    return(invisible(TRUE))
  }
  g0 <- tryCatch(as.numeric(gr(par)), error = function(e) NULL)
  if (is.null(g0) || length(g0) != length(par) || !all(is.finite(g0))) {
    return(invisible(TRUE))
  }
  nrm <- sqrt(sum(g0^2))
  if (nrm == 0) return(invisible(TRUE))
  d <- g0 / nrm
  h <- .Machine$double.eps^(1/3) * (1 + sqrt(sum(par^2)))
  fp <- tryCatch(fn(par + h * d), error = function(e) NA_real_)
  fm <- tryCatch(fn(par - h * d), error = function(e) NA_real_)
  if (length(fp) != 1L || length(fm) != 1L ||
      !is.finite(fp) || !is.finite(fm)) {
    return(invisible(TRUE))
  }
  # A stationary starting point carries no signal. The difference of `fn` along
  # a direction of no slope is its own truncation error, and comparing that
  # with a gradient of the same order is comparing two kinds of nothing: the
  # check fired on a caller who had done the best possible thing, handing in
  # the exact maximum likelihood estimate as the start. The scale a difference
  # can resolve is set by the size of the objective, not by 1.
  scale_f <- max(abs(fp), abs(fm), 1)
  if (nrm < 1e-6 * scale_f) return(invisible(TRUE))
  ana <- nrm                      # g'd with d = g/||g||
  off <- function(num) {
    abs(num - ana) / max(abs(num), abs(ana)) > 0.5 &&
      abs(num - ana) > 1e-6 * max(1, ana)
  }
  num <- (fp - fm) / (2 * h)
  # ONE STEP IS NOT A VERDICT. The step is scaled by the length of the whole
  # vector, so a coordinate in the hundreds sets it for every other one, and
  # near a stationary point a curved direction then reads its truncation
  # error as the slope: a score-driven filter on the Nile flows (an intercept
  # of 919 beside a log-loading) gave -0.00106 at h = 5.6e-3 against an exact
  # 0.00074, and agreed to four digits at h / 10. A wrong gradient disagrees
  # at every step, so the difference is repeated at two smaller steps before
  # anything is said.
  for (k in c(10, 100)) {
    if (!off(num)) break
    hk <- h / k
    fp <- tryCatch(fn(par + hk * d), error = function(e) NA_real_)
    fm <- tryCatch(fn(par - hk * d), error = function(e) NA_real_)
    if (length(fp) != 1L || length(fm) != 1L ||
        !is.finite(fp) || !is.finite(fm)) {
      break
    }
    num <- (fp - fm) / (2 * hk)
  }
  if (off(num)) {
    warning("'gr' does not appear to be the gradient of 'fn': along the ",
            "gradient direction at 'par',\n  'fn' changes at rate ",
            format(num, digits = 3), " where 'gr' predicts ",
            format(ana, digits = 3), ". Check that the two\n  compute the ",
            "same model. options(optimizers7.check_gradient = FALSE) turns ",
            "this check off.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Assemble the Result of a Run
#'
#' @description
#' Converts the list returned by the compiled loop into an
#' [optimizer_result()]: it fills in the point and the value, the counts and
#' the trace, sets the verdict, and composes the message.
#'
#' @details
#' `converged` is the flag returned by the loop, and it is never inferred
#' from the fact that the run ended. The loop sets it when the stopping rule
#' fires, when a line search with a declared resolution establishes that
#' every step fails to decrease the objective by more than that resolution,
#' and, for [brent()], when the bracket reaches the resolution of the
#' iterate. `criterion_met` is the label of the rule in the first case and a
#' label of its own in the others; a run that exhausted a budget reports the
#' budget.
#'
#' @param out The list from the C++ driver.
#' @param optimizer The optimizer that ran.
#' @param spec The objective handle, consulted for whether the gradient was
#'   supplied or differenced.
#' @param elapsed Seconds.
#' @param seed The generator state with which the run began, or `NULL` for a
#'   method that does not draw random numbers.
#'
#' @return An [optimizer_result()].
#'
#' @keywords internal
build_result <- function(out, optimizer, spec, elapsed, seed = NULL) {
  # Belt and braces beside the C++ fix: anything that is not an actual data
  # frame of rows becomes NULL, so that every consumer can test one thing.
  trace <- out$trace
  if (!is.data.frame(trace) || !nrow(trace)) trace <- NULL

  msg <- out$message
  # Only for a method that actually uses a gradient. A derivative-free one
  # differences nothing, and saying it did would be a false statement about how
  # exact the run was -- exactly the kind of claim this field exists to make
  # honestly. `gradient` comes back NULL from precisely those methods.
  if (!isTRUE(spec$has_gradient) && !is.null(out$gradient)) {
    extra <- "gradient obtained by finite differences"
    msg <- if (nzchar(msg)) paste(msg, extra, sep = "; ") else extra
  }

  stopped <- switch(out$stopped_by,
    criterion  = optimizer@criterion@label,
    maxit      = "iteration budget reached",
    max_eval   = "evaluation budget reached",
    failed     = "stopped without converging",
    resolution = "no decrease above the objective's resolution",
    bracket    = "bracket at machine resolution",
    out$stopped_by
  )

  optimizer_result(
    par = as.numeric(out$par),
    value = out$value,
    # NULL, not numeric(0): a derivative-free method computes no gradient, and
    # the class documents that absence as NULL. An empty numeric would print as
    # one and satisfy every is.null() test that exists to detect it.
    gradient = if (is.null(out$gradient)) NULL else as.numeric(out$gradient),
    counts = c(f = out$n_value, g = out$n_grad, h = out$n_hess),
    iterations = out$iterations,
    converged = isTRUE(out$converged),
    criterion_met = stopped,
    message = msg,
    trace = trace,
    optimizer = optimizer,
    elapsed = elapsed,
    seed = seed
  )
}


#' Run the Shared Descent Loop
#'
#' @description
#' The body of every method that has a direction: it prepares the objective,
#' passes the settings of the optimizer and the description of its direction
#' to the compiled loop, and assembles the result.
#'
#' @details
#' [newton()], [bfgs()], [lbfgs()], [gd()], [cg()], [bb()] and [chebyshev()]
#' differ only in the `method` list, which names the direction and carries
#' its parameters. The line search, the stopping rule, the budgets, the trace
#' and the reporting are the same code for all of them, so another method of
#' this kind needs a `Direction` in C++ and a constructor in R.
#'
#' @param optimizer The [optimizer()].
#' @param fn,par,gr,he The problem, as the user supplied it.
#' @param lower,upper Box constraints, as in [minimize()].
#' @param method A list describing the direction to the compiled loop.
#'
#' @return An [optimizer_result()].
#'
#' @keywords internal
run_descent <- function(optimizer, fn, par, gr, he, lower, upper, method) {
  spec <- prepare_objective(optimizer, fn, par, gr, he)
  bounds <- check_bounds(lower, upper, par)

  t0 <- proc.time()[["elapsed"]]
  out <- descent_run(
    spec = spec, par = as.numeric(par),
    criterion = optimizer@criterion, crit_fn = crit_met,
    method = method,
    line_search = line_search_spec(optimizer@line_search),
    maxit = as.integer(optimizer@maxit),
    max_eval = budget_int(optimizer@max_eval),
    verbose = optimizer@verbose,
    refresh = as.integer(optimizer@refresh),
    keep_trace = optimizer@keep_trace,
    step = optimizer@step,
    bounds = bounds
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  build_result(out, optimizer, spec, elapsed)
}


#' Validate an Initial Step Length
#'
#' @description
#' Checks that `step` is a single positive number, so that every
#' constructor that takes one rejects invalid values with the same message.
#'
#' @param step The value supplied.
#'
#' @return Invisibly `TRUE`. Raises an error naming `step` otherwise.
#'
#' @keywords internal
check_step <- function(step) {
  if (length(step) != 1L || !is.numeric(step) || is.na(step) || step <= 0) {
    stop("'step' must be a single positive number.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a Decrement Threshold
#'
#' @description
#' Checks that `stop_decrement` is a single positive number, or `Inf`, so that
#' [newton()] and [chebyshev()] reject invalid values with the same message.
#'
#' @param x The value supplied.
#'
#' @return Invisibly `TRUE`. Raises an error naming `stop_decrement` otherwise.
#'
#' @keywords internal
check_stop_decrement <- function(x) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || x <= 0) {
    stop("'stop_decrement' must be a single positive number (Inf to accept ",
         "every stop on a change).", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a Line Search
#'
#' @description
#' Checks that the value inherits from the abstract [line_search()] class.
#' The message names [armijo()] as an example, since passing the string
#' `"armijo"` instead of the object is a common mistake.
#'
#' @param x The value supplied.
#'
#' @return Invisibly `TRUE`. Raises an error naming `line_search` otherwise.
#'
#' @keywords internal
check_line_search <- function(x) {
  if (!S7::S7_inherits(x, line_search_class())) {
    stop("'line_search' must be a line_search object, e.g. armijo().",
         call. = FALSE)
  }
  invisible(TRUE)
}


#' @title Normalize Box Constraints
#'
#' @description
#' Recycles `lower` and `upper` to the length of `par`, checks
#' them, and returns them in the shape that the compiled loop reads.
#'
#' @details
#' The interface is two numeric vectors and not a list of pairs, as in
#' `stats::optim` and `stats::nlminb`, so that `lower = 0` makes every
#' parameter positive without one pair per coefficient. The list of pairs is
#' an internal shape, produced here, because the transform of each coordinate
#' in the compiled code reads it.
#'
#' Recycling is to length one or length `p` and nothing between: a `lower` of
#' length 2 for three parameters is much more likely to be a mistake than a
#' request, and the ordinary recycling of R would accept it silently.
#'
#' The starting value must be strictly interior. The reparametrization sends
#' a bound to an infinite value of the transformed variable, so a run that
#' starts exactly on a bound would begin at infinity, every later quantity
#' would be non-finite, and the failure would appear far from its cause. The
#' error raised here names the coordinate.
#'
#' @param lower,upper Numeric, of length one or `length(par)`.
#' @param par The starting value.
#'
#' @return A list with one `c(lower, upper)` pair per parameter, or an
#'   empty list when every bound is infinite. Call it at the top of a
#'   user-defined [minimize()] method, then hand each pair to
#'   [bounded_transform()]; an empty list means there is no box and
#'   the whole reparametrization should be skipped.
#'
#' @examples
#' check_bounds(0, Inf, par = c(1, 2))
#' check_bounds(c(0, -5), c(1, 5), par = c(0.5, 0))
#'
#' # infinite bounds on every parameter give no box at all
#' length(check_bounds(-Inf, Inf, par = c(1, 2)))
#'
#' # and a start on a bound is rejected, naming the coordinate
#' try(check_bounds(0, 1, par = c(0.5, 1)))
#'
#' @seealso [bounded_transform()], [minimize()]
#' @export
check_bounds <- function(lower, upper, par) {
  p <- length(par)
  fix <- function(v, nm) {
    if (is.null(v)) return(rep(if (nm == "lower") -Inf else Inf, p))
    if (!is.numeric(v) || anyNA(v)) {
      stop("'", nm, "' must be numeric and not NA.", call. = FALSE)
    }
    if (length(v) == 1L) return(rep(as.numeric(v), p))
    if (length(v) != p) {
      stop("'", nm, "' must have length 1 or ", p,
           ", one per parameter; it has length ", length(v), ".", call. = FALSE)
    }
    as.numeric(v)
  }
  lo <- fix(lower, "lower")
  up <- fix(upper, "upper")

  # Nothing finite means no box at all, and the compiled loop then skips the
  # whole reparametrization rather than composing with the identity p times.
  if (all(!is.finite(lo)) && all(!is.finite(up))) return(list())

  out <- vector("list", p)
  for (i in seq_len(p)) {
    if (lo[i] >= up[i]) {
      stop("For parameter ", i, " the lower bound must be strictly below the ",
           "upper one; they are ", format(lo[i]), " and ", format(up[i]), ".",
           call. = FALSE)
    }
    if (par[i] <= lo[i] || par[i] >= up[i]) {
      stop("The starting value for parameter ", i, " must lie strictly inside ",
           "its bounds (", format(lo[i]), ", ", format(up[i]), "); it is ",
           format(par[i]), ".", call. = FALSE)
    }
    out[[i]] <- c(lo[i], up[i])
  }
  out
}


#' @title The Bound Transform, and Its First Two Derivatives
#'
#' @description
#' Evaluates the reparametrization defined by a pair of bounds. The package
#' removes a box in this way, and the function is exported so that a
#' user-defined [minimize()] method can remove one in the same way.
#'
#' @param b A length-2 numeric vector, `c(lower, upper)`, using
#'   `-Inf` and `Inf` for a side that is unbounded.
#' @param eta A numeric vector on the unconstrained scale.
#'
#' @details
#' Bounds are removed here and not enforced. The map \eqn{h} is a shifted
#' exponential for a one-sided bound, a scaled logistic function for a
#' two-sided one, and the identity for neither (its inverse,
#' [bounded_forward()], is the shifted log or the scaled logit), so
#' optimizing in \eqn{\eta} makes every proposed point admissible by
#' construction.
#'
#' The image of \eqn{\eta = 0} is the middle of the box in the sense of the
#' map: `1` for \eqn{(0, \infty)}, `0.5` for \eqn{(0, 1)}, `4` for
#' \eqn{(2, 6)}, and `0` where there is no bound. For this reason
#' [start_zeros()] gives a sensible start for every kind of parameter.
#'
#' # Using it in a method written by the caller
#'
#' Map the starting value with [bounded_forward()], run unconstrained, and
#' wrap the objective so that it maps back before evaluating. The Jacobian is
#' diagonal, so the chain rule is one product per coordinate,
#' \eqn{\partial f/\partial \eta_i = (\partial f/\partial \theta_i)\, h_i'},
#' and `d2` is needed only when a Hessian is transformed, where it appears on
#' the diagonal alone. Report `par` on the caller's scale.
#'
#' # Relation to linkfunctions7
#'
#' \eqn{h} and its two derivatives are the inverse link of `bounded_link()`
#' of \pkg{linkfunctions7} and its derivatives, written in C++ because the transform is applied at every evaluation of the
#' objective, where a callback into R would cost more than the loop saves.
#' The test suite compares them with `linkinv()`, `dlinkinv()` and
#' `d2linkinv()`, so that the copy stays equal to the original.
#'
#' @return A named list of three numeric vectors, each as long as `eta`:
#'   \describe{
#'     \item{`h`}{the parameter \eqn{\theta = h(\eta)}, inside the box.}
#'     \item{`d1`}{\eqn{h'(\eta)}, the chain-rule factor for a gradient.}
#'     \item{`d2`}{\eqn{h''(\eta)}, needed only for a Hessian.}
#'   }
#'
#' @examples
#' # For a variance, the whole line maps onto the positive half-line, and every
#' # value of eta, however large, maps inside it.
#' bounded_transform(c(0, Inf), c(-2, 0, 2))$h
#' bounded_transform(c(0, Inf), c(-500, 500))$h > 0
#'
#' # For a probability, the result also includes the derivative that carries a
#' # gradient across.
#' str(bounded_transform(c(0, 1), c(-1, 0, 1)))
#'
#' # Zero maps to the middle of the box, whatever the box. For this reason
#' # start_zeros() gives a sensible start for every kind of parameter.
#' vapply(list(c(-Inf, Inf), c(0, Inf), c(0, 1), c(2, 6)),
#'        function(b) bounded_transform(b, 0)$h, 0)
#'
#' # The two maps are inverse to each other.
#' theta <- c(0.1, 0.5, 0.9)
#' all.equal(bounded_transform(c(0, 1), bounded_forward(c(0, 1), theta))$h,
#'           theta)
#'
#' # d1 is a derivative, and a central difference confirms it.
#' h <- 1e-6
#' (bounded_transform(c(0, 1), 0.3 + h)$h - bounded_transform(c(0, 1), 0.3 - h)$h) /
#'   (2 * h)
#' bounded_transform(c(0, 1), 0.3)$d1
#'
#' @seealso [bounded_forward()] for the inverse, [check_bounds()] for the
#'   shape in which `b` is supplied, [minimize()] for the run that uses both.
#' @export
bounded_transform <- function(b, eta) bounded_transform_cpp(b, as.numeric(eta))


#' @title The Forward Bound Transform
#'
#' @description
#' Maps a parameter from inside its bounds to the unconstrained scale: the
#' inverse of [bounded_transform()]'s `h`, and the map applied to a starting
#' value before an unconstrained run.
#'
#' @param b A length-2 numeric vector, `c(lower, upper)`.
#' @param theta A numeric vector strictly inside the bounds.
#'
#' @details
#' The values must lie strictly inside. A value on a bound maps to an
#' infinite \eqn{\eta}, so a run started there would begin at infinity and
#' fail far from its cause. For this reason [minimize()] rejects such a
#' starting value, with an error that names the coordinate, before anything
#' is evaluated.
#'
#' @return A numeric vector on the unconstrained scale, as long as `theta`.
#'   A value on or outside a bound gives an infinite or `NaN` entry rather
#'   than an error; the check belongs to [check_bounds()].
#'
#' @examples
#' bounded_forward(c(0, Inf), c(0.5, 1, 8))
#' bounded_forward(c(0, 1), 0.5)
#'
#' # A value on a bound has no finite image, so such a starting point is
#' # unusable and minimize() rejects it.
#' bounded_forward(c(0, 1), c(0, 1))
#'
#' # The round trip also holds in this direction.
#' eta <- c(-2, 0, 3)
#' all.equal(bounded_forward(c(0, Inf), bounded_transform(c(0, Inf), eta)$h),
#'           eta)
#'
#' @seealso [bounded_transform()] for the inverse and its derivatives,
#'   [check_bounds()] for the validation.
#' @export
bounded_forward <- function(b, theta) bounded_forward_cpp(b, as.numeric(theta))
