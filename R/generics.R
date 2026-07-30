#' @include optimizer_class.R
#' @include result_class.R
#' @include objective.R
NULL

#' @title Minimise a Function
#'
#' @description
#' The entry point of the package. Everything here minimises; see
#' \code{\link{maximize}} for the other direction.
#'
#' @param optimizer An \code{\link{optimizer}} object, carrying the algorithm and
#'   its settings.
#' @param fn The objective: an ordinary R function of the parameter vector, a
#'   \code{\link{finite_sum}} object, or a \code{\link{cpp_objective}}.
#' @param par A numeric vector of starting values.
#' @param gr An optional gradient function. Ignored when \code{fn} carries its
#'   own gradient.
#' @param he An optional Hessian function. Methods that do not use one ignore
#'   it, so calling code need not branch on the algorithm.
#' @param ... Passed to methods.
#'
#' @details
#' The generic dispatches on \code{optimizer} alone, so each algorithm is
#' written once. Which of the three shapes \code{fn} arrived in is settled
#' separately by \code{\link{as_objective}}, which dispatches on \code{fn} and
#' hands every algorithm the same handle.
#'
#' A gradient that is not supplied is computed by central finite differences.
#' The result records which derivatives were supplied and which were
#' differenced, so a run is never silently less exact than it appears.
#'
#' @return An \code{\link{optimizer_result}}.
#'
#' @examples
#' # a quadratic, with the gradient supplied
#' minimize(gradient_descent(), function(p) sum((p - c(1, 2))^2),
#'          par = c(0, 0), gr = function(p) 2 * (p - c(1, 2)))
#'
#' # and without, so the gradient is differenced
#' minimize(gradient_descent(), function(p) sum((p - c(1, 2))^2), c(0, 0))
#'
#' @seealso \code{\link{maximize}}, \code{\link{gradient_descent}},
#'   \code{\link{crit_any}}
#' @export
minimize <- S7::new_generic("minimize", "optimizer",
  function(optimizer, fn, par, gr = NULL, he = NULL, ...) S7::S7_dispatch())


#' @title Maximise a Function
#'
#' @description
#' Runs \code{\link{minimize}} on the negated objective and reports the value
#' with its sign restored.
#'
#' @inheritParams minimize
#'
#' @details
#' Every algorithm in the package minimises, always; that is the convention and
#' the names say so. This is the thin wrapper for the other direction, and it
#' exists so that nobody has to remember to negate a log-likelihood by hand and
#' then negate the answer back.
#'
#' @return An \code{\link{optimizer_result}}, whose \code{value} and
#'   \code{gradient} refer to the original objective rather than the negated one.
#'
#' @examples
#' maximize(gradient_descent(), function(p) -sum((p - c(1, 2))^2), c(0, 0))
#'
#' @seealso \code{\link{minimize}}
#' @export
maximize <- function(optimizer, fn, par, gr = NULL, he = NULL, ...) {
  if (!is.function(fn)) {
    stop("maximize() takes a plain function; negate other objectives yourself.",
         call. = FALSE)
  }
  neg_fn <- function(p) -fn(p)
  neg_gr <- if (is.null(gr)) NULL else function(p) -gr(p)
  neg_he <- if (is.null(he)) NULL else function(p) -he(p)

  res <- minimize(optimizer, neg_fn, par, gr = neg_gr, he = neg_he, ...)

  res@value <- -res@value
  if (!is.null(res@gradient)) res@gradient <- -res@gradient
  if (!is.null(res@trace)) res@trace$value <- -res@trace$value
  res
}
