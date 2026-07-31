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


#' @title Normalise an Objective for the Optimisers
#'
#' @description
#' Turns whatever the user supplied — an R function or a compiled C++ function
#' behind an external pointer — into the single handle the algorithms are
#' written against.
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
#' @seealso \code{\link{cpp_objective}}
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


# --- compiled objectives ----------------------------------------------------

#' @title A Compiled C++ Objective
#'
#' @description
#' Marks a pair of external pointers as an objective the optimisers may call
#' directly, so that the iteration never returns to R.
#'
#' @param fn An external pointer to a C++ function
#'   \code{double f(const arma::vec&)}.
#' @param gr An optional external pointer to
#'   \code{arma::vec g(const arma::vec&)}.
#'
#' @details
#' This is what makes writing the algorithms in C++ worth the trouble. With an R
#' objective every evaluation is a callback into R, and on a typical problem the
#' objective dominates so thoroughly that the compiled loop saves little. With a
#' compiled objective the loop never leaves C++.
#'
#' The pointers are \strong{tagged} with an S3 class rather than passed bare.
#' A bare external pointer may be anything at all — a database handle, a pointer
#' to an unrelated object — and calling through the wrong one is not an R error
#' but a crash of the process. Requiring a tag means the package only calls
#' pointers that something deliberately marked as callable.
#'
#' @return An object of class \code{cpp_objective}.
#'
#' @examples
#' # The pointers come from compiled code, so the example shows the guard
#' # rather than the happy path: anything that is not an external pointer is
#' # refused, because calling through the wrong one crashes the process
#' # instead of raising an error.
#' try(cpp_objective(function(p) sum(p^2)))
#'
#' @seealso \code{\link{as_objective}}, \code{\link{minimize}}
#' @export
cpp_objective <- function(fn, gr = NULL) {
  # typeof() is the whole test: an external pointer has no class attribute to
  # inherit from, so inherits() and methods::is() add nothing here except a
  # dependency and a way to be wrong.
  if (typeof(fn) != "externalptr") {
    stop("'fn' must be an external pointer to a C++ function.", call. = FALSE)
  }
  if (!is.null(gr) && typeof(gr) != "externalptr") {
    stop("'gr' must be an external pointer or NULL.", call. = FALSE)
  }
  structure(list(fn = fn, gr = gr), class = "cpp_objective")
}


#' @title A Compiled Objective, as an Objective
#' @name as_objective.cpp_objective
#' @description
#' Passes the pointers through to the C++ side, which calls them without
#' returning to R.
#' @param fn A \code{\link{cpp_objective}}.
#' @param gr Ignored; the gradient belongs to the object.
#' @param ... Unused.
#' @return An objective handle; see \code{\link{as_objective}}.
#' @keywords internal
S7::method(as_objective, S7::new_S3_class("cpp_objective")) <-
  function(fn, gr = NULL, he = NULL, ...) {
    list(kind = "cpp", fn_ptr = fn$fn, gr_ptr = fn$gr,
         has_gradient = !is.null(fn$gr), has_hessian = FALSE)
  }
