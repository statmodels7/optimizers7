#' @include optimizer_class.R
#' @include generics.R
NULL

#' @title S7 Class for Gradient Descent
#' @description The class \code{\link{gradient_descent}} instantiates.
#' @param step The initial step length offered to the line search.
#' @param line_search A \code{\link{line_search}} object.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{gradient_descent}}
#' @keywords internal
GradientDescent <- S7::new_class("GradientDescent", parent = optimizer,
  properties = list(
    step        = S7::class_numeric,
    line_search = S7::class_any
  ))


#' @title Gradient Descent with Backtracking
#'
#' @description
#' The simplest descent method: step along the negative gradient, shrinking the
#' step until the objective improves.
#'
#' @param criterion The stopping rule; see \code{\link{crit_any}}.
#' @param step Initial step length offered to the line search each iteration.
#'   Defaults to \code{1}.
#' @param line_search How far to go along the direction; see
#'   \code{\link{armijo}} and \code{\link{wolfe}}. Defaults to Armijo
#'   backtracking, which is all a first-order method needs.
#' @param maxit Maximum iterations. Defaults to 500.
#' @param max_eval Maximum objective evaluations. Defaults to 10000.
#' @param verbose Report progress? Defaults to \code{FALSE}.
#' @param refresh Report every this many iterations. Defaults to 10.
#' @param keep_trace Store the iteration path? Defaults to \code{FALSE}.
#'
#' @details
#' Nobody should reach for this to solve a real problem — it converges linearly
#' and slowly, and \code{bfgs()} will beat it on anything smooth. It is here
#' because it is the smallest algorithm that exercises the whole frame: the
#' objective abstraction in all three of its shapes, the stopping rule evaluated
#' as an object, the trace, the reporting, the evaluation budget, and the
#' refusal to step to a point where the objective is not finite.
#'
#' That last one is the safeguard worth naming. A step is accepted only if it
#' both improves the objective and produces a finite value; otherwise it is
#' shrunk and tried again. Propagating one \code{NaN} contaminates every iterate
#' after it, and the run then fails somewhere far from the cause.
#'
#' @return An S7 object of class \code{GradientDescent}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @examples
#' gradient_descent()
#' gradient_descent(criterion = crit_grad(1e-10), maxit = 2000)
#'
#' @seealso \code{\link{minimize}}, \code{\link{crit_any}}
#' @export
gradient_descent <- function(criterion = crit_any(crit_grad(1e-8),
                                                  crit_rel_obj(1e-12)),
                             step = 1, line_search = armijo(),
                             maxit = 500, max_eval = 10000,
                             verbose = FALSE, refresh = 10,
                             keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  if (length(step) != 1L || !is.numeric(step) || is.na(step) || step <= 0) {
    stop("'step' must be a single positive number.", call. = FALSE)
  }
  if (!S7::S7_inherits(line_search, line_search_class())) {
    stop("'line_search' must be a line_search object, e.g. armijo().",
         call. = FALSE)
  }
  GradientDescent(
    name = "gradient descent", criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    step = step, line_search = line_search
  )
}


#' @title Minimise by Gradient Descent
#' @name minimize.GradientDescent
#' @description
#' Runs \code{\link{gradient_descent}} on the objective.
#' @param optimizer A \code{GradientDescent} object.
#' @param fn,par,gr,he,bounds,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, GradientDescent) <-
  function(optimizer, fn, par, gr = NULL, he = NULL, bounds = NULL, ...) {
    run_descent(optimizer, fn, par, gr, he, bounds,
                list(type = "gradient_descent"))
  }
