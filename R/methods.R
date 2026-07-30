#' @include generics.R
#' @include result_class.R
NULL

# The two pieces every algorithm's method needs, written once here so that a new
# algorithm is a C++ loop and a constructor and nothing else.


#' Prepare an Objective for an Algorithm
#'
#' @description
#' Normalises the objective, checks the stopping rule can actually be evaluated
#' by this optimiser, and checks the starting value.
#'
#' @details
#' The criterion check is the interesting one. A rule needing a gradient handed
#' to a method that computes none would sit there testing \code{NULL} at every
#' iteration and never fire, so the run would end on the iteration budget and
#' report failure for a reason nowhere near the truth. Refusing it here, by
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

  needs <- crit_needs(optimizer@criterion)
  provides <- optimizer_provides(optimizer)
  missing <- setdiff(needs, provides)
  if (length(missing)) {
    stop(sprintf(
      paste0("The stopping rule needs %s, which %s does not compute.\n",
             "  Choose a criterion this optimiser can evaluate, or a method ",
             "that provides it."),
      paste(missing, collapse = ", "), optimizer@name), call. = FALSE)
  }
  spec
}


#' What an Optimiser Can Offer a Stopping Rule
#'
#' @description
#' The names of the \code{state} components an optimiser is able to fill in, so
#' that \code{\link{prepare_objective}} can refuse a rule it could never satisfy.
#'
#' @details
#' Gradient-based methods provide a gradient; the derivative-free ones to come
#' will not, and will override this to say so.
#'
#' @param optimizer An \code{\link{optimizer}}.
#'
#' @return A character vector.
#'
#' @keywords internal
optimizer_provides <- S7::new_generic("optimizer_provides", "optimizer",
                                      function(optimizer) S7::S7_dispatch())

S7::method(optimizer_provides, optimizer) <- function(optimizer) "gradient"


#' Assemble the Result of a Run
#'
#' @description
#' Turns what the C++ loop returned into an \code{\link{optimizer_result}}.
#'
#' @details
#' The one judgement here is the meaning of \code{converged}: it is taken
#' straight from whether the stopping rule fired, and never inferred from the
#' run having ended. An optimiser that exhausted its iterations has not
#' converged, and saying otherwise turns a failure into a wrong answer that
#' looks right.
#'
#' @param out The list from the C++ driver.
#' @param optimizer The optimiser that ran.
#' @param spec The objective handle, consulted for whether the gradient was
#'   supplied or differenced.
#' @param elapsed Seconds.
#'
#' @return An \code{\link{optimizer_result}}.
#'
#' @keywords internal
build_result <- function(out, optimizer, spec, elapsed) {
  # Belt and braces beside the C++ fix: anything that is not an actual data
  # frame of rows becomes NULL, so that every consumer can test one thing.
  trace <- out$trace
  if (!is.data.frame(trace) || !nrow(trace)) trace <- NULL

  msg <- out$message
  if (!isTRUE(spec$has_gradient)) {
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
    gradient = as.numeric(out$gradient),
    counts = c(f = out$n_value, g = out$n_grad, h = out$n_hess),
    iterations = out$iterations,
    converged = isTRUE(out$converged),
    criterion_met = stopped,
    message = msg,
    trace = trace,
    optimizer = optimizer,
    elapsed = elapsed
  )
}


#' Run the Shared Descent Loop
#'
#' @description
#' The body of every method that has a direction: prepares the objective, hands
#' the compiled loop the optimiser's settings and the description of its
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
#' @param method A list describing the direction to the compiled loop.
#'
#' @return An \code{\link{optimizer_result}}.
#'
#' @keywords internal
run_descent <- function(optimizer, fn, par, gr, he, method) {
  spec <- prepare_objective(optimizer, fn, par, gr, he)

  t0 <- proc.time()[["elapsed"]]
  out <- descent_run(
    spec = spec, par = as.numeric(par),
    criterion = optimizer@criterion, crit_fn = crit_met,
    method = method,
    line_search = line_search_spec(optimizer@line_search),
    maxit = as.integer(optimizer@maxit),
    max_eval = as.integer(optimizer@max_eval),
    verbose = optimizer@verbose,
    refresh = as.integer(optimizer@refresh),
    keep_trace = optimizer@keep_trace,
    step = optimizer@step
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
