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
#' @param note An optional extra sentence appended when a criterion is refused,
#'   for a method whose reason is not simply that it computes no such quantity.
#'
#' @return The objective handle from \code{\link{as_objective}}.
#'
#' @keywords internal
prepare_objective <- function(optimizer, fn, par, gr = NULL, he = NULL,
                              note = NULL) {
  if (!is.numeric(par) || !length(par) || anyNA(par)) {
    stop("'par' must be a numeric vector of starting values.", call. = FALSE)
  }
  spec <- as_objective(fn, gr, he)

  needs <- crit_needs(optimizer@criterion)
  provides <- optimizer_provides(optimizer)
  missing <- setdiff(needs, provides)
  if (length(missing)) {
    msg <- sprintf(
      paste0("The stopping rule needs %s, which %s does not provide.\n",
             "  Choose a criterion this optimiser can evaluate, or a method ",
             "that provides it."),
      paste(missing, collapse = ", "), optimizer@name)
    if (!is.null(note)) msg <- paste0(msg, "\n  ", note)
    stop(msg, call. = FALSE)
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
#' Gradient-based methods provide a gradient and a reproducible objective value;
#' the derivative-free ones to come will not provide the first, and a stochastic
#' method provides neither, since both are then estimates whose noise a
#' tolerance would be measuring instead of the progress.
#'
#' @param optimizer An \code{\link{optimizer}}.
#'
#' @return A character vector.
#'
#' @keywords internal
optimizer_provides <- S7::new_generic("optimizer_provides", "optimizer",
                                      function(optimizer) S7::S7_dispatch())

S7::method(optimizer_provides, optimizer) <- function(optimizer)
  c("gradient", "objective")


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
#' @param bounds Optional box constraints, as in \code{\link{minimize}}.
#' @param method A list describing the direction to the compiled loop.
#'
#' @return An \code{\link{optimizer_result}}.
#'
#' @keywords internal
run_descent <- function(optimizer, fn, par, gr, he, bounds, method) {
  spec <- prepare_objective(optimizer, fn, par, gr, he)
  bounds <- check_bounds(bounds, par)

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


#' Validate Box Constraints
#'
#' @description
#' Checks that the bounds are well formed and that the starting value lies
#' strictly inside them, and returns them in the shape the compiled loop reads.
#'
#' @details
#' The starting value must be strictly interior, and refusing a boundary start
#' is not pedantry. The reparametrisation sends a bound to an infinite value of
#' the transformed variable, so a run started exactly on one begins at infinity:
#' every subsequent quantity is non-finite and the failure surfaces far from its
#' cause. Saying so here, naming the coordinate, costs one check.
#'
#' @param bounds A list of length-2 numeric vectors, or \code{NULL}.
#' @param par The starting value.
#'
#' @return A list of length-2 numeric vectors, or an empty list when there are
#'   no bounds — the compiled side reads an empty list as "unconstrained".
#'
#' @keywords internal
check_bounds <- function(bounds, par) {
  if (is.null(bounds)) return(list())
  if (!is.list(bounds) || length(bounds) != length(par)) {
    stop("'bounds' must be a list with one element per parameter.",
         call. = FALSE)
  }
  out <- vector("list", length(par))
  for (i in seq_along(bounds)) {
    b <- bounds[[i]]
    if (!is.numeric(b) || length(b) != 2L || anyNA(b)) {
      stop("Element ", i, " of 'bounds' must be c(lower, upper).", call. = FALSE)
    }
    if (b[1] >= b[2]) {
      stop("In 'bounds' element ", i, ", the lower bound must be strictly ",
           "below the upper one.", call. = FALSE)
    }
    if (par[i] <= b[1] || par[i] >= b[2]) {
      stop("The starting value for parameter ", i, " must lie strictly inside ",
           "its bounds (", format(b[1]), ", ", format(b[2]), "); it is ",
           format(par[i]), ".", call. = FALSE)
    }
    out[[i]] <- as.numeric(b)
  }
  out
}


#' The Bound Transform, for Checking Against linkfunctions7
#'
#' @description
#' Evaluates the reparametrisation a set of bounds implies, together with its
#' first two derivatives.
#'
#' @details
#' Not part of the interface. It exists because the transforms are written out
#' in C++ rather than called through \pkg{linkfunctions7} — the transform is
#' applied on every objective evaluation, and a callback into R there would undo
#' the reason for compiling the loop — so something has to hold the copy against
#' the original. The test suite does, on every run.
#'
#' @param b A length-2 numeric vector, \code{c(lower, upper)}.
#' @param eta A numeric vector on the unconstrained scale.
#'
#' @return A list with \code{h}, \code{d1} and \code{d2}, matching
#'   \code{linkinv()}, \code{dlinkinv()} and \code{d2linkinv()} of the
#'   corresponding \code{bounded_link()}.
#'
#' @keywords internal
bounded_transform <- function(b, eta) bounded_transform_cpp(b, as.numeric(eta))


#' The Forward Bound Transform
#'
#' @description
#' Maps a parameter inside its bounds to the unconstrained scale; the inverse of
#' \code{\link{bounded_transform}}'s \code{h}.
#'
#' @param b A length-2 numeric vector, \code{c(lower, upper)}.
#' @param theta A numeric vector strictly inside the bounds.
#'
#' @return A numeric vector on the unconstrained scale.
#'
#' @keywords internal
bounded_forward <- function(b, theta) bounded_forward_cpp(b, as.numeric(theta))
