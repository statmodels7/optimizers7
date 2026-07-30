#' @include criterion.R
NULL

# The objective, and the three shapes it may arrive in.
#
# minimize() dispatches on the optimizer, so an algorithm is written once.
# as_objective() dispatches on the objective, so the three shapes are told apart
# once. Each generic dispatches where there is real variation; dispatching
# minimize() on both would need one method per algorithm per shape, and every
# algorithm would be written three times.


#' @title Normalise an Objective for the Optimisers
#'
#' @description
#' Turns whatever the user supplied — an R function, a finite sum over
#' observations, or a compiled C++ function behind an external pointer — into the
#' single handle the algorithms are written against.
#'
#' @param fn The objective.
#' @param gr An optional gradient.
#' @param ... Passed to methods.
#'
#' @return A list describing the objective to the C++ side: \code{kind}, the
#'   pieces belonging to that kind, and \code{has_gradient}.
#'
#' @seealso \code{\link{finite_sum}}, \code{\link{cpp_objective}}
#' @export
as_objective <- S7::new_generic("as_objective", "fn",
                                function(fn, gr = NULL, ...) S7::S7_dispatch())


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
S7::method(as_objective, S7::class_function) <- function(fn, gr = NULL, ...) {
  if (!is.null(gr) && !is.function(gr)) {
    stop("'gr' must be a function or NULL.", call. = FALSE)
  }
  list(kind = "r", fn = fn, gr = gr, has_gradient = !is.null(gr))
}


# --- finite sums ------------------------------------------------------------

#' @title S7 Class for a Finite-Sum Objective
#'
#' @description
#' The class \code{\link{finite_sum}} instantiates: an objective that is a sum
#' over observations and can therefore be evaluated on a subset of them.
#'
#' @param fn A function \code{fn(par, idx)}.
#' @param gr An optional \code{gr(par, idx)}.
#' @param n The number of terms.
#'
#' @return An S7 object of class \code{FiniteSum}.
#'
#' @seealso \code{\link{finite_sum}}
#' @keywords internal
FiniteSum <- S7::new_class("FiniteSum",
  properties = list(
    fn = S7::class_function,
    gr = S7::class_any,
    n  = S7::class_numeric
  ))


#' @title An Objective That Is a Sum Over Observations
#'
#' @description
#' Declares that the objective is \eqn{\sum_{i=1}^{n} f_i(\theta)}, so that a
#' stochastic method may evaluate it on a subset of the terms.
#'
#' @param fn A function \code{fn(par, idx)} returning the sum of the terms
#'   indexed by \code{idx}.
#' @param gr An optional \code{gr(par, idx)} returning that sum's gradient.
#' @param n The number of terms.
#'
#' @details
#' Subsampling is not a property of an algorithm; it is a property of the
#' objective. A black-box \code{fn(par)} has no terms to draw from — the
#' optimiser does not know what an observation is and must not guess — so a
#' method that wants to subsample can only be given an objective that has said
#' it may. That is what this constructor is for.
#'
#' The full-sample value is \code{fn(par, seq_len(n))}, so an ordinary
#' evaluation is a special case and nothing is duplicated.
#'
#' @return An S7 object of class \code{FiniteSum}, accepted by
#'   \code{\link{minimize}} anywhere an objective is.
#'
#' @examples
#' set.seed(1)
#' y <- rnorm(100, mean = 2)
#'
#' # the negative Gaussian log-likelihood in the mean, as a sum over observations
#' obj <- finite_sum(
#'   fn = function(par, idx) sum((y[idx] - par)^2) / 2,
#'   gr = function(par, idx) -sum(y[idx] - par),
#'   n  = length(y)
#' )
#' obj@n
#'
#' @seealso \code{\link{as_objective}}, \code{\link{minimize}}
#' @export
finite_sum <- function(fn, gr = NULL, n) {
  if (!is.function(fn)) stop("'fn' must be a function.", call. = FALSE)
  if (!is.null(gr) && !is.function(gr)) {
    stop("'gr' must be a function or NULL.", call. = FALSE)
  }
  if (length(n) != 1L || !is.numeric(n) || is.na(n) || n < 1) {
    stop("'n' must be a single positive number of terms.", call. = FALSE)
  }
  if (length(formals(fn)) < 2L) {
    stop("'fn' must take two arguments, fn(par, idx).", call. = FALSE)
  }
  FiniteSum(fn = fn, gr = gr, n = as.integer(n))
}


#' @title A Finite Sum as an Objective
#' @name as_objective.FiniteSum
#' @description
#' Wraps the two-argument functions so that the algorithms see the usual
#' one-argument interface, and keeps the pieces a stochastic method needs to
#' subsample.
#' @param fn A \code{\link{FiniteSum}} object.
#' @param gr Ignored; the gradient belongs to the object.
#' @param ... Unused.
#' @return An objective handle; see \code{\link{as_objective}}.
#' @keywords internal
S7::method(as_objective, FiniteSum) <- function(fn, gr = NULL, ...) {
  obj <- fn
  n <- obj@n
  all_idx <- seq_len(n)
  gr_obj <- obj@gr
  list(
    kind = "finite_sum",
    fn = function(par) obj@fn(par, all_idx),
    gr = if (is.null(gr_obj)) NULL else function(par) gr_obj(par, all_idx),
    fn_idx = obj@fn,
    gr_idx = gr_obj,
    n = n,
    has_gradient = !is.null(gr_obj)
  )
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
  function(fn, gr = NULL, ...) {
    list(kind = "cpp", fn_ptr = fn$fn, gr_ptr = fn$gr,
         has_gradient = !is.null(fn$gr))
  }
