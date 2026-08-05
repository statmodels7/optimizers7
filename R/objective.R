#' @include criterion.R
NULL

# The objective, and the two shapes it may arrive in.
#
# minimize() dispatches on the optimizer, so an algorithm is written once.
# as_objective() dispatches on the objective, so the two shapes are told apart
# once. Each generic dispatches where there is real variation; dispatching
# minimize() on both would need one method per algorithm per shape, and every
# algorithm would be written twice.
#
# There was a third, finite_sum(), an objective declaring itself a sum over
# observations so that a stochastic method could ask it for a subsample. It
# existed for exactly one caller, Adam, and went when Adam stopped drawing its
# own minibatches: an objective that resamples is a closure, and needs no class.


#' @title Normalize an Objective for the Optimizers
#'
#' @description
#' Turns the objective into the single handle the algorithms are written
#' against.
#'
#' @param fn The objective.
#' @param gr An optional gradient.
#' @param he An optional Hessian. Only Newton uses one; the other methods
#'   accept it and ignore it, so calling code need not branch on the method.
#' @param ... Passed to methods.
#'
#' @return A list describing the objective to the C++ side: \code{kind}, the
#'   pieces belonging to that kind, and flags saying which derivatives were
#'   supplied rather than differenced.
#'
#' @examples
#' # a plain function, with and without a gradient
#' str(as_objective(function(p) sum(p^2)))
#' str(as_objective(function(p) sum(p^2), gr = function(p) 2 * p))
#'
#' @export
as_objective <- S7::new_generic("as_objective", "fn",
                                function(fn, gr = NULL, he = NULL, ...)
                                  S7::S7_dispatch())


#' @title An Ordinary R Function as an Objective
#' @name as_objective.function
#' @description
#' The common case: \code{fn(par)} returns a number and, if supplied,
#' \code{gr(par)} returns the gradient.
#' @param fn A function of the parameter vector.
#' @param gr An optional gradient function, or \code{NULL} for finite differences.
#' @param ... Unused.
#' @return An objective handle; see \code{\link{as_objective}}.
#' @keywords internal
S7::method(as_objective, S7::class_function) <-
  function(fn, gr = NULL, he = NULL, ...) {
    if (!is.null(gr) && !is.function(gr)) {
      stop("'gr' must be a function or NULL.", call. = FALSE)
    }
    if (!is.null(he) && !is.function(he)) {
      stop("'he' must be a function or NULL.", call. = FALSE)
    }
    list(kind = "r", fn = fn, gr = gr, he = he,
         has_gradient = !is.null(gr), has_hessian = !is.null(he))
  }
