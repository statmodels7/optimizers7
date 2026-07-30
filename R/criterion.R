#' @title S7 Class for Convergence Criteria
#'
#' @import S7
#' @description
#' A stopping rule, as an object. Every optimiser carries one, and the user may
#' replace it, combine several, or write a new kind.
#'
#' @details
#' The alternative would be an argument taking a string, and a `switch` inside
#' every algorithm. That is exactly the arrangement this toolkit exists to
#' replace: it fixes the set of rules at the moment the package is written, and
#' nothing outside can add to it. A criterion here is an object implementing one
#' generic, \code{\link{crit_met}}, so a rule the authors never thought of is a
#' first-class citizen.
#'
#' Criteria are combined with \code{\link{crit_any}} and \code{\link{crit_all}},
#' which are themselves criteria, so combinations nest.
#'
#' @param label A short character label, used when reporting which rule fired.
#'
#' @return An S7 object of class \code{criterion}. The class is abstract: use one
#'   of the constructors, or write your own subclass.
#'
#' @seealso \code{\link{crit_grad}}, \code{\link{crit_rel_obj}},
#'   \code{\link{crit_any}}, \code{\link{crit_met}}
#' @export
criterion <- S7::new_class(
  "criterion",
  properties = list(label = S7::class_character),
  abstract = TRUE
)


#' @title Has the Stopping Rule Been Met?
#'
#' @description
#' The one generic a criterion must implement.
#'
#' @param criterion A \code{\link{criterion}} object.
#' @param state A named list describing the current iteration; see Details.
#'
#' @details
#' \code{state} carries everything any rule could need:
#' \describe{
#'   \item{\code{iter}}{the iteration just completed.}
#'   \item{\code{f_new}, \code{f_old}}{the objective after and before it.}
#'   \item{\code{x_new}, \code{x_old}}{the parameter vectors, likewise.}
#'   \item{\code{gradient}}{the gradient at \code{x_new}, or \code{NULL} when the
#'     method does not compute one.}
#' }
#' A rule that needs something absent from \code{state} — a gradient, from a
#' derivative-free method — must say so through \code{\link{crit_needs}} rather
#' than silently never firing.
#'
#' @return A single logical.
#'
#' @examples
#' st <- list(iter = 3, f_new = 1.0000001, f_old = 1.0000002,
#'            x_new = c(1, 2), x_old = c(1, 2), gradient = c(1e-9, -2e-9))
#' crit_met(crit_grad(1e-8), st)
#' crit_met(crit_abs_obj(1e-12), st)
#'
#' @export
crit_met <- S7::new_generic("crit_met", "criterion",
                            function(criterion, state) S7::S7_dispatch())


#' @title What a Criterion Needs From the Iteration
#'
#' @description
#' The names of the \code{state} components a criterion requires, so that an
#' algorithm can refuse a rule it cannot evaluate instead of accepting one that
#' never fires.
#'
#' @details
#' A derivative-free method has no gradient, so \code{\link{crit_grad}} would sit
#' there testing \code{NULL} at every iteration and quietly never stop the run.
#' Refusing it at construction is the same discipline that makes
#' \code{check_link()} in \pkg{linkfunctions7} report a numerical derivative
#' order as numerical rather than as passed: a check that cannot be evaluated
#' must say so.
#'
#' @param criterion A \code{\link{criterion}} object.
#'
#' @return A character vector, possibly empty.
#'
#' @examples
#' crit_needs(crit_grad(1e-8))
#' crit_needs(crit_rel_obj(1e-10))
#'
#' @export
crit_needs <- S7::new_generic("crit_needs", "criterion",
                              function(criterion) S7::S7_dispatch())

S7::method(crit_needs, criterion) <- function(criterion) character()


# --- gradient ---------------------------------------------------------------

#' @title S7 Class for the Gradient Criterion
#' @description The class \code{\link{crit_grad}} instantiates.
#' @param tol The tolerance.
#' @param norm Either \code{"max"} or \code{"2"}.
#' @return An S7 object inheriting from \code{\link{criterion}}.
#' @seealso \code{\link{crit_grad}}
#' @keywords internal
CritGrad <- S7::new_class("CritGrad", parent = criterion,
  properties = list(tol = S7::class_numeric, norm = S7::class_character))

S7::method(crit_needs, CritGrad) <- function(criterion) "gradient"

S7::method(crit_met, CritGrad) <- function(criterion, state) {
  g <- state$gradient
  if (is.null(g) || !length(g) || anyNA(g)) return(FALSE)
  v <- if (identical(criterion@norm, "2")) sqrt(sum(g^2)) else max(abs(g))
  v < criterion@tol
}

#' @title Stop When the Gradient Is Small
#'
#' @description
#' The rule \eqn{\lVert \nabla f \rVert < \texttt{tol}}.
#'
#' @param tol Numeric tolerance. Defaults to \code{1e-8}.
#' @param norm \code{"max"} (default) or \code{"2"}.
#'
#' @details
#' The max-norm is the default because it does not grow with the dimension the
#' way the 2-norm does: the same tolerance then means the same thing for a
#' two-parameter problem and a two-hundred-parameter one, whereas \code{1e-8} in
#' the 2-norm is a far stricter demand in high dimension. Both are available, so
#' the choice is only a default.
#'
#' Only usable by a method that computes a gradient; a derivative-free optimiser
#' refuses it rather than accepting a rule that can never fire.
#'
#' @return A \code{\link{criterion}} object.
#'
#' @examples
#' crit_grad()
#' crit_grad(1e-10, norm = "2")
#'
#' @seealso \code{\link{crit_rel_obj}}, \code{\link{crit_any}}
#' @export
crit_grad <- function(tol = 1e-8, norm = c("max", "2")) {
  norm <- match.arg(norm)
  check_tol(tol)
  CritGrad(label = paste0("gradient (", norm, "-norm) < ", format(tol)),
           tol = tol, norm = norm)
}


# --- objective --------------------------------------------------------------

#' @title S7 Class for the Absolute Objective Criterion
#' @description The class \code{\link{crit_abs_obj}} instantiates.
#' @param tol The tolerance.
#' @return An S7 object inheriting from \code{\link{criterion}}.
#' @seealso \code{\link{crit_abs_obj}}
#' @keywords internal
CritAbsObj <- S7::new_class("CritAbsObj", parent = criterion,
  properties = list(tol = S7::class_numeric))

S7::method(crit_met, CritAbsObj) <- function(criterion, state) {
  if (is.null(state$f_old) || !is.finite(state$f_old)) return(FALSE)
  abs(state$f_new - state$f_old) < criterion@tol
}

#' @title Stop When the Objective Stops Moving (Absolute)
#'
#' @description
#' The rule \eqn{\lvert f_{new} - f_{old} \rvert < \texttt{tol}}.
#'
#' @param tol Numeric tolerance. Defaults to \code{1e-10}.
#'
#' @return A \code{\link{criterion}} object.
#'
#' @examples
#' crit_abs_obj()
#'
#' @seealso \code{\link{crit_rel_obj}}
#' @export
crit_abs_obj <- function(tol = 1e-10) {
  check_tol(tol)
  CritAbsObj(label = paste0("|df| < ", format(tol)), tol = tol)
}


#' @title S7 Class for the Relative Objective Criterion
#' @description The class \code{\link{crit_rel_obj}} instantiates.
#' @param tol The tolerance.
#' @return An S7 object inheriting from \code{\link{criterion}}.
#' @seealso \code{\link{crit_rel_obj}}
#' @keywords internal
CritRelObj <- S7::new_class("CritRelObj", parent = criterion,
  properties = list(tol = S7::class_numeric))

S7::method(crit_met, CritRelObj) <- function(criterion, state) {
  if (is.null(state$f_old) || !is.finite(state$f_old)) return(FALSE)
  # The floor is not decoration: without it an optimum sitting at zero divides
  # by zero and the rule either never fires or fires at once.
  abs(state$f_new - state$f_old) <
    criterion@tol * (abs(state$f_old) + criterion@tol)
}

#' @title Stop When the Objective Stops Moving (Relative)
#'
#' @description
#' The rule
#' \eqn{\lvert f_{new} - f_{old} \rvert < \texttt{tol}\,(\lvert f_{old} \rvert + \texttt{tol})}.
#'
#' @param tol Numeric tolerance. Defaults to \code{1e-12}.
#'
#' @details
#' The \code{+ tol} in the denominator is a floor, and it is load-bearing: an
#' objective whose optimum is at zero would otherwise be compared against a
#' vanishing scale.
#'
#' @return A \code{\link{criterion}} object.
#'
#' @examples
#' crit_rel_obj()
#'
#' @seealso \code{\link{crit_abs_obj}}
#' @export
crit_rel_obj <- function(tol = 1e-12) {
  check_tol(tol)
  CritRelObj(label = paste0("|df| < ", format(tol), " (relative)"), tol = tol)
}


# --- parameters -------------------------------------------------------------

#' @title S7 Class for the Absolute Parameter Criterion
#' @description The class \code{\link{crit_abs_par}} instantiates.
#' @param tol The tolerance.
#' @return An S7 object inheriting from \code{\link{criterion}}.
#' @seealso \code{\link{crit_abs_par}}
#' @keywords internal
CritAbsPar <- S7::new_class("CritAbsPar", parent = criterion,
  properties = list(tol = S7::class_numeric))

S7::method(crit_met, CritAbsPar) <- function(criterion, state) {
  if (is.null(state$x_old)) return(FALSE)
  max(abs(state$x_new - state$x_old)) < criterion@tol
}

#' @title Stop When the Parameters Stop Moving (Absolute)
#'
#' @description
#' The rule \eqn{\max_j \lvert x_j^{new} - x_j^{old} \rvert < \texttt{tol}}.
#'
#' @param tol Numeric tolerance. Defaults to \code{1e-8}.
#'
#' @return A \code{\link{criterion}} object.
#'
#' @examples
#' crit_abs_par()
#'
#' @seealso \code{\link{crit_rel_par}}
#' @export
crit_abs_par <- function(tol = 1e-8) {
  check_tol(tol)
  CritAbsPar(label = paste0("|dx| < ", format(tol)), tol = tol)
}


#' @title S7 Class for the Relative Parameter Criterion
#' @description The class \code{\link{crit_rel_par}} instantiates.
#' @param tol The tolerance.
#' @return An S7 object inheriting from \code{\link{criterion}}.
#' @seealso \code{\link{crit_rel_par}}
#' @keywords internal
CritRelPar <- S7::new_class("CritRelPar", parent = criterion,
  properties = list(tol = S7::class_numeric))

S7::method(crit_met, CritRelPar) <- function(criterion, state) {
  if (is.null(state$x_old)) return(FALSE)
  max(abs(state$x_new - state$x_old) /
        (abs(state$x_old) + criterion@tol)) < criterion@tol
}

#' @title Stop When the Parameters Stop Moving (Relative)
#'
#' @description
#' The rule
#' \eqn{\max_j \lvert x_j^{new} - x_j^{old}\rvert / (\lvert x_j^{old}\rvert + \texttt{tol}) < \texttt{tol}}.
#'
#' @param tol Numeric tolerance. Defaults to \code{1e-8}.
#'
#' @return A \code{\link{criterion}} object.
#'
#' @examples
#' crit_rel_par()
#'
#' @seealso \code{\link{crit_abs_par}}
#' @export
crit_rel_par <- function(tol = 1e-8) {
  check_tol(tol)
  CritRelPar(label = paste0("|dx| < ", format(tol), " (relative)"), tol = tol)
}


# --- combinators ------------------------------------------------------------

#' @title S7 Class for a Combination of Criteria
#' @description The class \code{\link{crit_any}} and \code{\link{crit_all}}
#'   instantiate.
#' @param criteria A list of \code{\link{criterion}} objects.
#' @param how Either \code{"any"} or \code{"all"}.
#' @return An S7 object inheriting from \code{\link{criterion}}.
#' @seealso \code{\link{crit_any}}, \code{\link{crit_all}}
#' @keywords internal
CritCombine <- S7::new_class("CritCombine", parent = criterion,
  properties = list(criteria = S7::class_list, how = S7::class_character))

S7::method(crit_needs, CritCombine) <- function(criterion) {
  unique(unlist(lapply(criterion@criteria, crit_needs)))
}

S7::method(crit_met, CritCombine) <- function(criterion, state) {
  met <- vapply(criterion@criteria, crit_met, logical(1), state = state)
  if (identical(criterion@how, "all")) all(met) else any(met)
}

# Shared body of the two combinators.
combine_criteria <- function(dots, how) {
  if (!length(dots)) {
    stop("At least one criterion is required.", call. = FALSE)
  }
  ok <- vapply(dots, function(z) S7::S7_inherits(z, criterion), logical(1))
  if (!all(ok)) {
    stop("Every argument must be a 'criterion' object.", call. = FALSE)
  }
  sep <- if (identical(how, "all")) " and " else " or "
  CritCombine(
    label = paste(vapply(dots, function(z) z@label, character(1)), collapse = sep),
    criteria = dots, how = how
  )
}

#' @title Stop When Any of Several Rules Fires
#'
#' @description
#' Combines criteria disjunctively. This is the usual arrangement: a run should
#' end as soon as any reasonable rule is satisfied.
#'
#' @param ... \code{\link{criterion}} objects.
#'
#' @return A \code{\link{criterion}} object, so combinations nest.
#'
#' @examples
#' crit_any(crit_grad(1e-8), crit_rel_obj(1e-12))
#'
#' @seealso \code{\link{crit_all}}
#' @export
crit_any <- function(...) combine_criteria(list(...), "any")

#' @title Stop Only When Every Rule Fires
#'
#' @description
#' Combines criteria conjunctively, for a run that should not stop until several
#' independent things agree.
#'
#' @param ... \code{\link{criterion}} objects.
#'
#' @return A \code{\link{criterion}} object, so combinations nest.
#'
#' @examples
#' crit_all(crit_grad(1e-6), crit_abs_par(1e-10))
#'
#' @seealso \code{\link{crit_any}}
#' @export
crit_all <- function(...) combine_criteria(list(...), "all")


#' @title Print Method for Criteria
#' @name print.criterion
#' @param x A \code{\link{criterion}} object.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @examples
#' print(crit_any(crit_grad(), crit_rel_obj()))
#' @keywords internal
S7::method(print, criterion) <- function(x, ...) {
  cat("<criterion> ", x@label, "\n", sep = "")
  invisible(x)
}


#' Validate a Tolerance
#'
#' @description
#' Every criterion constructor takes a tolerance, and every one of them should
#' reject the same nonsense in the same words.
#'
#' @param tol The value supplied.
#'
#' @return Invisibly \code{TRUE}; raises an error otherwise.
#'
#' @keywords internal
check_tol <- function(tol) {
  if (length(tol) != 1L || !is.numeric(tol) || is.na(tol) || tol <= 0) {
    stop("'tol' must be a single positive number.", call. = FALSE)
  }
  invisible(TRUE)
}
