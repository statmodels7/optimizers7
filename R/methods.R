#' @include generics.R
#' @include result_class.R
NULL

# The two pieces every algorithm's method needs, written once here so that a new
# algorithm is a C++ loop and a constructor and nothing else.


#' Prepare an Objective for an Algorithm
#'
#' @description
#' Normalizes the objective, checks the stopping rule can actually be evaluated
#' by this optimizer, and checks the starting value.
#'
#' @details
#' The criterion check is the interesting one. A rule needing a gradient handed
#' to a method that computes none would sit there testing \code{NULL} at every
#' iteration and never fire, so the run would end on the iteration budget and
#' report failure for a reason nowhere near the truth. Rejecting it here, by
#' name, is the same discipline as \code{check_link()} in \pkg{linkfunctions7}
#' reporting a numerical derivative order as numerical rather than as passed: a
#' check that cannot be evaluated must say so rather than pass or fail silently.
#'
#' @param optimizer The \code{\link{optimizer}}.
#' @param fn,gr,he The objective and its optional derivatives, as supplied.
#' @param par The starting value.
#'
#' @return The objective handle from \code{\link{as_objective}}.
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


#' @title Reject a Stopping Rule an Optimizer Cannot Evaluate
#'
#' @description
#' Compares what the optimizer's criterion reads against what the optimizer can
#' supply, and stops if the rule asks for something absent.
#'
#' @param optimizer An \code{\link{optimizer}}.
#'
#' @details
#' Call this at the top of a user-defined \code{\link{minimize}} method. A rule
#' needing a gradient, handed to a method that computes none, would sit testing
#' \code{NULL} at every iteration and never fire; the run would then end on its
#' iteration budget and report a reason nowhere near the truth. Rejecting it here,
#' by name, is what \code{\link{check_optimizer}} tests for.
#'
#' What an optimizer can supply is declared by
#' \code{\link{optimizer_provides}}, whose default claims a gradient. Override
#' that when the method computes none.
#'
#' @return Invisibly \code{TRUE}; raises an error otherwise.
#'
#' @examples
#' check_criterion(bfgs())
#' try(check_criterion(nelder_mead(criterion = crit_grad())))
#'
#' @seealso \code{\link{optimizer_provides}}, \code{\link{crit_needs}}
#' @export
check_criterion <- function(optimizer) {
  needs <- crit_needs(optimizer@criterion)
  provides <- optimizer_provides(optimizer)
  missing <- setdiff(needs, provides)
  if (length(missing)) {
    msg <- sprintf(
      paste0("The stopping rule needs %s, which %s does not provide.\n",
             "  Choose a criterion this optimizer can evaluate, or a method ",
             "that provides it."),
      paste(missing, collapse = ", "), optimizer@name)
    stop(msg, call. = FALSE)
  }
  invisible(TRUE)
}


#' @title What an Optimizer Can Offer a Stopping Rule
#'
#' @description
#' The names of the \code{state} components an optimizer is able to fill in, so
#' that \code{\link{check_criterion}} can reject a rule it could never satisfy.
#'
#' @details
#' Gradient-based methods provide a gradient. The derivative-free ones provide
#' a stationarity measure instead, since no single derivative they could report
#' goes to zero at a solution.
#'
#' Every optimizer evaluates the objective, so there is no token for that and
#' rules reading it are never rejected. A user-defined method inherits the
#' default, which claims a gradient; if that is not true of it, say so, because
#' the rejection machinery relies on this declaration being accurate.
#'
#' @param optimizer An \code{\link{optimizer}}.
#'
#' @return A character vector. The names the shipped criteria read are
#'   \code{"gradient"} and \code{"stationarity"}.
#'
#' @examples
#' optimizer_provides(bfgs())
#' optimizer_provides(nelder_mead())
#'
#' @seealso \code{\link{check_criterion}}, \code{\link{crit_needs}}
#' @export
optimizer_provides <- S7::new_generic("optimizer_provides", "optimizer",
                                      function(optimizer) S7::S7_dispatch())

S7::method(optimizer_provides, optimizer) <- function(optimizer) "gradient"


#' @title Whether an Optimizer Takes Box Bounds
#'
#' @description
#' \code{TRUE} when the optimizer honours the \code{lower} and \code{upper}
#' arguments of \code{\link{minimize}}, so that
#' \code{\link{check_optimizer}} tests them, and \code{FALSE} for a method
#' that takes its constraint another way.
#'
#' @details
#' Every method of the package removes bounds by reparametrization and
#' answers \code{TRUE}, which is the default. A proximal method is the
#' exception: a constraint reaches it inside the proximal operator, where it
#' composes with the term already there, so bounds beside the objective would
#' be a second and conflicting route to the same thing.
#'
#' @param optimizer An \code{\link{optimizer}}.
#'
#' @return A single logical.
#'
#' @examples
#' optimizer_bounded(bfgs())
#'
#' @seealso \code{\link{optimizer_provides}}
#' @export
optimizer_bounded <- S7::new_generic("optimizer_bounded", "optimizer",
                                     function(optimizer) S7::S7_dispatch())

S7::method(optimizer_bounded, optimizer) <- function(optimizer) TRUE


#' Translate an Evaluation Budget for the Compiled Loop
#'
#' @description
#' The compiled loops hold the budget as an integer, and the default budget is
#' \code{Inf}: no cap at all. This maps \code{Inf} to the largest representable
#' integer, which no run reaches, and any finite value to itself.
#'
#' @param x The \code{max_eval} property, a single positive number or \code{Inf}.
#' @return A single integer.
#' @keywords internal
budget_int <- function(x) {
  if (is.finite(x)) as.integer(x) else .Machine$integer.max
}


#' Warn When a Supplied Gradient Does Not Belong to the Objective
#'
#' @description
#' Compares the directional derivative of \code{fn} along the gradient
#' direction at \code{par}, obtained from one central difference, with the rate
#' the supplied \code{gr} predicts. A gradient computed from a different model
#' than the objective -- a fixed predictor where the objective uses the
#' parameter, a stale copy of the data -- makes the two disagree grossly, and
#' the symptom downstream is otherwise a mute line-search failure at the first
#' iteration.
#'
#' @details
#' The check costs one call to \code{gr} and two to \code{fn}, runs once per
#' \code{\link{minimize}} call, and is skipped whenever it cannot be decisive:
#' a non-function objective, a non-finite gradient at \code{par}, an objective
#' that is not finite at the probe points, and a gradient too small to be
#' distinguished from the truncation error of the difference it is compared
#' with --- which is what a caller who starts at the optimum supplies. The
#' tolerance is
#' deliberately loose -- a relative disagreement above one half -- so that
#' finite-difference error or a subgradient of a non-smooth objective does not
#' trip it; it exists to catch the wrong function, not the eighth digit.
#'
#' Setting \code{options(optimizers7.check_gradient = FALSE)} disables it, and
#' \code{\link{multistart}} disables it for its inner runs after the first.
#'
#' @param fn The objective.
#' @param gr The supplied gradient.
#' @param par The starting point, numeric.
#'
#' @return Invisibly \code{TRUE}; warns on a gross mismatch.
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
  num <- (fp - fm) / (2 * h)
  ana <- nrm                      # g'd with d = g/||g||
  rel <- abs(num - ana) / max(abs(num), abs(ana))
  if (rel > 0.5 && abs(num - ana) > 1e-6 * max(1, ana)) {
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
#' Turns what the C++ loop returned into an \code{\link{optimizer_result}}.
#'
#' @details
#' The one judgement here is the meaning of \code{converged}: it is taken
#' straight from whether the stopping rule fired, and never inferred from the
#' run having ended. An optimizer that exhausted its iterations has not
#' converged, and saying otherwise turns a failure into a wrong answer that
#' looks right.
#'
#' @param out The list from the C++ driver.
#' @param optimizer The optimizer that ran.
#' @param spec The objective handle, consulted for whether the gradient was
#'   supplied or differenced.
#' @param elapsed Seconds.
#' @param seed The generator state the run began with, or \code{NULL} for a
#'   method that draws no random numbers.
#'
#' @return An \code{\link{optimizer_result}}.
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
    criterion = optimizer@criterion@label,
    maxit     = "iteration budget reached",
    max_eval  = "evaluation budget reached",
    failed    = "stopped without converging",
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
#' The body of every method that has a direction: prepares the objective, hands
#' the compiled loop the optimizer's settings and the description of its
#' direction, and assembles the result.
#'
#' @details
#' Newton, BFGS, L-BFGS and gradient descent differ only in the \code{method}
#' list, which names the direction and carries its parameters. Everything
#' else — the line search, the stopping rule, the budgets, the trace, the
#' reporting — is the same code for all of them, which is what makes adding a
#' fifth method a Direction in C++ and a constructor in R.
#'
#' @param optimizer The \code{\link{optimizer}}.
#' @param fn,par,gr,he The problem, as the user supplied it.
#' @param lower,upper Box constraints, as in \code{\link{minimize}}.
#' @param method A list describing the direction to the compiled loop.
#'
#' @return An \code{\link{optimizer_result}}.
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
#' @param step The value supplied.
#' @return Invisibly \code{TRUE}; raises an error otherwise.
#' @keywords internal
check_step <- function(step) {
  if (length(step) != 1L || !is.numeric(step) || is.na(step) || step <= 0) {
    stop("'step' must be a single positive number.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a Line Search
#' @param x The value supplied.
#' @return Invisibly \code{TRUE}; raises an error otherwise.
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
#' Recycles \code{lower} and \code{upper} to the length of \code{par}, checks
#' them, and returns them in the shape the compiled loop reads.
#'
#' @details
#' The interface is two numeric vectors rather than a list of pairs, which is
#' what \code{stats::optim} and \code{stats::nlminb} take and what makes
#' \code{lower = 0} say "every parameter is positive" in four characters instead
#' of a list of identical pairs one per coefficient. The list of pairs is an
#' internal shape, produced here, because that is what the per-coordinate
#' transform on the C++ side wants.
#'
#' Recycling is length one or length \code{p} and nothing between: a
#' \code{lower} of length 2 for three parameters is far more likely to be a
#' mistake than a request, and R's ordinary recycling would silently oblige.
#'
#' The starting value must be strictly interior, and rejecting a boundary start
#' is not pedantry. The reparametrization sends a bound to an infinite value of
#' the transformed variable, so a run started exactly on one begins at infinity:
#' every subsequent quantity is non-finite and the failure surfaces far from its
#' cause. Saying so here, naming the coordinate, costs one check.
#'
#' @param lower,upper Numeric, of length one or \code{length(par)}.
#' @param par The starting value.
#'
#' @return A list with one \code{c(lower, upper)} pair per parameter, or an
#'   empty list when no bound is finite. Call it at the top of a
#'   user-defined \code{\link{minimize}} method, then hand each pair to
#'   \code{\link{bounded_transform}}; an empty list means there is no box and
#'   the whole reparametrization should be skipped.
#'
#' @examples
#' check_bounds(0, Inf, par = c(1, 2))
#' check_bounds(c(0, -5), c(1, 5), par = c(0.5, 0))
#'
#' # no finite bound is no box at all
#' length(check_bounds(-Inf, Inf, par = c(1, 2)))
#'
#' # and a start on a bound is rejected, naming the coordinate
#' try(check_bounds(0, 1, par = c(0.5, 1)))
#'
#' @seealso \code{\link{bounded_transform}}, \code{\link{minimize}}
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
#' Evaluates the reparametrization a set of bounds implies. This is how the
#' package removes a box, and it is exported so that a \code{\link{minimize}}
#' user-defined method can remove one the same way.
#'
#' @param b A length-2 numeric vector, \code{c(lower, upper)}, using
#'   \code{-Inf} and \code{Inf} for a side that is unbounded.
#' @param eta A numeric vector on the unconstrained scale.
#'
#' @details
#' Bounds are not enforced here, they are removed: a shifted log for a one-sided
#' bound, a scaled logit for two, and the identity for neither. Optimize in
#' \eqn{\eta} and every proposed point is admissible by construction.
#'
#' To honor \code{bounds} in a user-defined method: map the starting value with
#' \code{\link{bounded_forward}}, run unconstrained, and wrap the objective so
#' that it maps back before evaluating. The Jacobian is diagonal, so the chain
#' rule is short —
#' \eqn{\partial f/\partial \eta_i = (\partial f/\partial \theta_i)\, h_i'} —
#' and \code{d2} is needed only when a Hessian is transformed, where it appears on
#' the diagonal alone. Report \code{par} on the user's scale.
#'
#' These are \pkg{linkfunctions7}'s \code{bounded_link()}, written out in C++
#' because the transform is applied on every objective evaluation and a callback
#' into R there would undo the reason for compiling the loop. The test suite
#' pins them to \code{linkinv()}, \code{dlinkinv()} and \code{d2linkinv()} on
#' every run, so the copy cannot drift from the original.
#'
#' @return A list with \code{h} (the parameter), \code{d1} and \code{d2}.
#'
#' @examples
#' # a variance: the whole line maps onto the positive half
#' bounded_transform(c(0, Inf), c(-2, 0, 2))$h
#'
#' # a probability, and the derivative that carries a gradient across
#' str(bounded_transform(c(0, 1), c(-1, 0, 1)))
#'
#' # the round trip
#' bounded_transform(c(0, 1), bounded_forward(c(0, 1), c(0.1, 0.5, 0.9)))$h
#'
#' @seealso \code{\link{bounded_forward}}, \code{\link{minimize}}
#' @export
bounded_transform <- function(b, eta) bounded_transform_cpp(b, as.numeric(eta))


#' @title The Forward Bound Transform
#'
#' @description
#' Maps a parameter from inside its bounds to the unconstrained scale: the
#' inverse of \code{\link{bounded_transform}}'s \code{h}, and what a starting
#' value goes through before an unconstrained run.
#'
#' @param b A length-2 numeric vector, \code{c(lower, upper)}.
#' @param theta A numeric vector strictly inside the bounds.
#'
#' @details
#' Strictly inside. A value on a bound maps to an infinite \eqn{\eta}, so a run
#' started there begins at infinity and fails far from its cause; that is why
#' \code{\link{minimize}} rejects such a starting value by name.
#'
#' @return A numeric vector on the unconstrained scale.
#'
#' @examples
#' bounded_forward(c(0, Inf), c(0.5, 1, 8))
#' bounded_forward(c(0, 1), 0.5)
#'
#' @seealso \code{\link{bounded_transform}}
#' @export
bounded_forward <- function(b, theta) bounded_forward_cpp(b, as.numeric(theta))
