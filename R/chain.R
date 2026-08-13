#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
#' @include check_optimizer.R
NULL

#' @title S7 Class for a Sequence of Optimizers
#' @description The class \code{\link{chain}} instantiates.
#' @param stages The optimizers, in the order they run.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{chain}}
#' @name Chain-class
#' @aliases Chain
#' @keywords internal
Chain <- S7::new_class("Chain", parent = optimizer,
  properties = list(stages = S7::class_list))


#' @title Run Optimizers One After Another
#'
#' @description
#' Wraps several optimizers into one, each starting where the previous
#' finished.
#'
#' @details
#' The composition a global search needs: \code{chain(sa(), lbfgs())} explores
#' first and then descends from wherever the exploration left off, and neither
#' method has to know about the other. It is the second wrapper of this shape
#' after \code{\link{multistart}}, and the two compose --
#' \code{multistart(chain(sa(), lbfgs()))} is a legal optimizer.
#'
#' Each stage carries its OWN criterion and its own budgets, which is the point
#' of a chain rather than an argument: a coarse rule and a small budget for the
#' exploration, a tight rule for the descent.
#'
#' \strong{What the result reports.} The point and the value are the LAST
#' stage's, since that is where the run ended, and so is \code{converged}: a
#' chain has converged when the method that finished it says so, and an earlier
#' stage exhausting its own budget is the ordinary way a global search ends
#' rather than a failure of the whole. Evaluations and iterations are summed
#' over the stages, and the trace, when kept, carries a \code{stage} column.
#'
#' A stage that raises propagates, since a method that cannot run on the
#' objective is a fact about the objective. A stage that runs without
#' converging passes its point on, which is what the first stage of a chain
#' usually does.
#'
#' \code{chain(x)} with a single stage is that stage's run, reported through
#' the chain: the wrapper is not a special case to be avoided.
#'
#' @param ... Two or more optimizers, in the order they should run. A single
#'   one is accepted.
#' @param verbose Report which stage is running? Defaults to \code{FALSE}.
#' @param keep_trace Store the path of every stage? Defaults to \code{FALSE}.
#'
#' @return An S7 object of class \code{Chain}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @examples
#' chain(sa(maxit = 10), bfgs())
#'
#' # a multimodal objective: the search says which basin, the descent finishes
#' # the job inside it
#' rastrigin <- function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p))
#' set.seed(3)
#' minimize(sa(maxit = 20), rastrigin, c(3.5, -2.5))@value
#' set.seed(3)
#' minimize(chain(sa(maxit = 20), bfgs()), rastrigin, c(3.5, -2.5))@value
#'
#' @seealso \code{\link{multistart}}, \code{\link{sa}}
#' @export
chain <- function(..., verbose = FALSE, keep_trace = FALSE) {
  stages <- list(...)
  if (!length(stages)) {
    stop("'chain' needs at least one optimizer.", call. = FALSE)
  }
  ok <- vapply(stages, function(o) S7::S7_inherits(o, optimizer_class()),
               logical(1))
  if (!all(ok)) {
    stop("every stage must be an optimizer object, e.g. sa() or bfgs().",
         call. = FALSE)
  }
  if (length(verbose) != 1L || !is.logical(verbose) || is.na(verbose)) {
    stop("'verbose' must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(keep_trace) != 1L || !is.logical(keep_trace) ||
      is.na(keep_trace)) {
    stop("'keep_trace' must be TRUE or FALSE.", call. = FALSE)
  }
  last <- stages[[length(stages)]]
  Chain(
    name = paste(vapply(stages, function(o) o@name, character(1)),
                 collapse = " then "),
    # The rule and the budgets shown are the last stage's, because the result
    # is reported from there: check_criterion() then asks the question that
    # decides whether the reported run could evaluate its own rule.
    criterion = last@criterion, maxit = last@maxit, max_eval = last@max_eval,
    verbose = verbose, refresh = 1, keep_trace = keep_trace,
    stages = stages
  )
}


#' @title What a Chain Can Offer a Stopping Rule
#' @name optimizer_provides.Chain
#' @description
#' Whatever the LAST stage offers, that being the stage whose rule ends the run
#' and whose result is reported.
#' @param optimizer A \code{Chain} object.
#' @return A character vector.
#' @keywords internal
S7::method(optimizer_provides, Chain) <- function(optimizer)
  optimizer_provides(optimizer@stages[[length(optimizer@stages)]])


#' @title Whether a Chain Takes Box Bounds
#' @name optimizer_bounded.Chain
#' @description
#' \code{TRUE} only when EVERY stage takes them: bounds are passed to all of
#' them, so one stage that would ignore them makes the chain unable to promise
#' the box.
#' @param optimizer A \code{Chain} object.
#' @return A single logical.
#' @keywords internal
S7::method(optimizer_bounded, Chain) <- function(optimizer)
  all(vapply(optimizer@stages, optimizer_bounded, logical(1)))


#' @rdname with_criterion
#' @name with_criterion.Chain
#' @keywords internal
S7::method(with_criterion, Chain) <- function(optimizer, criterion) {
  # The rule that decides the run belongs to the last stage; setting only the
  # outer copy would change the printing and nothing else.
  st <- optimizer@stages
  st[[length(st)]] <- with_criterion(st[[length(st)]], criterion)
  S7::set_props(optimizer, criterion = criterion, stages = st)
}

#' @rdname with_maxit
#' @name with_maxit.Chain
#' @keywords internal
S7::method(with_maxit, Chain) <- function(optimizer, maxit) {
  st <- optimizer@stages
  st[[length(st)]] <- with_maxit(st[[length(st)]], maxit)
  S7::set_props(optimizer, maxit = maxit, stages = st)
}


#' @title Minimize by a Sequence of Optimizers
#' @name minimize.Chain
#' @description Runs each stage from where the previous one finished.
#' @param optimizer A \code{Chain} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}, passed to
#'   every stage.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Chain) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    stages <- optimizer@stages
    # Every stage's rule is checked before any of them runs, so a chain whose
    # third stage cannot evaluate its criterion fails before spending the first
    # two rather than after.
    for (st in stages) check_criterion(st)

    # The generic checked fn against gr once already, at the caller's par;
    # letting each stage repeat it would print the same warning once per stage.
    old_opt <- options(optimizers7.check_gradient = FALSE)
    on.exit(options(old_opt), add = TRUE)

    x <- as.numeric(par)
    t0 <- proc.time()[["elapsed"]]
    results <- vector("list", length(stages))
    for (k in seq_along(stages)) {
      if (optimizer@verbose) {
        cat(sprintf("  stage %d/%d: %s\n", k, length(stages),
                    stages[[k]]@name))
      }
      res <- minimize(stages[[k]], fn, x, gr = gr, he = he,
                      lower = lower, upper = upper, ...)
      results[[k]] <- res
      x <- res@par
    }
    elapsed <- proc.time()[["elapsed"]] - t0

    build_chain_result(results, optimizer, elapsed)
  }


#' Assemble the Result of a Chain
#'
#' @description
#' The last stage's point and verdict, with the work of every stage summed and
#' the traces stacked.
#'
#' @param results The per-stage results, in order.
#' @param optimizer The \code{Chain}.
#' @param elapsed Total seconds.
#'
#' @return An \code{\link{optimizer_result}}.
#'
#' @keywords internal
build_chain_result <- function(results, optimizer, elapsed) {
  last <- results[[length(results)]]
  # counts is a named c(f =, g =, h =): the work of the whole chain is what a
  # caller comparing it against a single optimizer wants to see. Added rather
  # than summed coordinate by coordinate so that the STORAGE TYPE survives --
  # a chain of one stage must report that stage's counts identically, and
  # sum() over a vapply would hand back a double where the stage had integers.
  counts <- Reduce(`+`, lapply(results, function(r) r@counts))

  tr <- NULL
  if (optimizer@keep_trace) {
    parts <- list()
    for (k in seq_along(results)) {
      t_k <- results[[k]]@trace
      if (is.data.frame(t_k) && nrow(t_k)) {
        t_k$stage <- k
        parts[[length(parts) + 1L]] <- t_k
      }
    }
    # the stages need not report the same columns: a derivative-free one has a
    # stationarity where a descent has a gradient norm, so they are stacked
    # only when they agree and otherwise the last stage's is kept
    if (length(parts)) {
      same <- length(unique(lapply(parts, names))) == 1L
      tr <- if (same) do.call(rbind, parts) else parts[[length(parts)]]
    }
  }

  msg <- vapply(seq_along(results), function(k) {
    m <- results[[k]]@message
    if (nzchar(m)) sprintf("stage %d: %s", k, m) else ""
  }, character(1))
  msg <- paste(msg[nzchar(msg)], collapse = "; ")

  optimizer_result(
    par = last@par, value = last@value, gradient = last@gradient,
    counts = counts, iterations = last@iterations,
    converged = last@converged, criterion_met = last@criterion_met,
    message = msg, trace = tr, optimizer = optimizer, elapsed = elapsed,
    seed = results[[1L]]@seed
  )
}
