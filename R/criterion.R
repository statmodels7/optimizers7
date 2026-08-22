#' @title S7 Class for Convergence Criteria
#'
#' @import S7
#' @description
#' A stopping rule, as an object. Every optimizer carries one, and the user may
#' replace it, combine several, or write a new kind.
#'
#' @details
#' The alternative would be an argument taking a string, and a `switch` inside
#' every algorithm. That is exactly the arrangement this toolkit exists to
#' replace: it fixes the set of rules at the moment the package is written, and
#' nothing outside can add to it. A criterion here is an object implementing one
#' generic, \code{\link{crit_met}}, so a user-defined rule is treated like
#' a shipped one.
#'
#' Criteria are combined with \code{\link{crit_any}} and \code{\link{crit_all}},
#' which are themselves criteria, so combinations nest.
#'
#' @param label A short character label, used when reporting which rule fired.
#'
#' @return An S7 object of class \code{criterion}. The class is abstract: use one
#'   of the constructors, or from a user-defined subclass.
#'
#' @examples
#' # The class is abstract, so it cannot be instantiated directly...
#' try(criterion(label = "mine"))
#'
#' # ...but anything inheriting from it is a criterion, including a rule the
#' # package never anticipated.
#' Tiny <- S7::new_class("Tiny", parent = criterion,
#'                       properties = list(tol = S7::class_numeric))
#' S7::method(crit_met, Tiny) <- function(criterion, state)
#'   state$f_new < criterion@tol
#' crit_met(Tiny(label = "f < 1e-6", tol = 1e-6), list(f_new = 1e-9))
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
#'   \item{\code{stationarity}}{a non-negative measure of remaining progress,
#'     supplied by the derivative-free methods in place of a gradient, or
#'     \code{NULL}. See \code{\link{crit_stationary}}.}
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
#' @seealso \code{\link{as_objective}}, \code{\link{crit_needs}}, \code{\link{check_criterion}}
#' @export
crit_met <- S7::new_generic("crit_met", "criterion",
                            function(criterion, state) S7::S7_dispatch())


#' @title What a Criterion Needs From the Iteration
#'
#' @description
#' The names of the \code{state} components a criterion requires, so that an
#' algorithm can reject a rule it cannot evaluate instead of accepting one that
#' never fires.
#'
#' @details
#' A derivative-free method has no gradient, so \code{\link{crit_grad}} would sit
#' there testing \code{NULL} at every iteration and quietly never stop the run.
#' Rejecting it at construction is the same discipline that makes
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
#' @seealso \code{\link{as_objective}}, \code{\link{crit_met}}, \code{\link{check_criterion}}
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
#' @param tol Numeric tolerance. Defaults to \code{1e-6}.
#' @param norm \code{"max"} (default) or \code{"2"}.
#'
#' @details
#' The max-norm is the default because it does not grow with the dimension the
#' way the 2-norm does: the same tolerance then means the same thing for a
#' two-parameter problem and a two-hundred-parameter one, whereas \code{1e-6} in
#' the 2-norm is a far stricter demand in high dimension. Both are available, so
#' the choice is only a default.
#'
#' How small a gradient a run can actually reach is set by the objective, not
#' by the method. A line search accepts a step only when the objective
#' decreases by a definite amount, and near a minimum that decrease is about
#' \eqn{\lVert \nabla f \rVert^{2} / (2\lambda)} for a curvature \eqn{\lambda}.
#' Once it drops below the rounding of the objective itself, about
#' \eqn{\varepsilon \lvert f \rvert}, no step in any direction can be verified
#' and the search stops, so the smallest attainable gradient is around
#' \eqn{\sqrt{2 \lambda \varepsilon \lvert f^{*} \rvert}} and grows with the
#' value at the solution. On conjugate gradients applied to Rosenbrock, adding
#' a constant to the objective --- which moves neither the minimizer nor the
#' gradient --- takes the attainable gradient from \code{1.9e-9} at
#' \eqn{f^{*} = 0} to \code{4.4e-8} at \eqn{f^{*} = 1} and \code{6.5e-5} at
#' \eqn{f^{*} = 10^{6}}. The default suits an objective of order one at its
#' solution, which is what a log-likelihood per observation is; an objective
#' that lands in the millions needs a correspondingly looser tolerance, and one
#' that lands at zero can be asked for much more.
#'
#' Only usable by a method that computes a gradient; a derivative-free optimizer
#' rejects it rather than accepting a rule that can never fire.
#'
#' @return A \code{\link{criterion}} object.
#'
#' @examples
#' crit_grad()
#' crit_grad(1e-10, norm = "2")
#'
#' @seealso \code{\link{crit_rel_obj}}, \code{\link{crit_any}}
#' @export
crit_grad <- function(tol = 1e-6, norm = c("max", "2")) {
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


# --- stationarity -----------------------------------------------------------

#' @title S7 Class for the Stationarity Criterion
#' @description The class \code{\link{crit_stationary}} instantiates.
#' @param tol The tolerance.
#' @return An S7 object inheriting from \code{\link{criterion}}.
#' @seealso \code{\link{crit_stationary}}
#' @keywords internal
CritStationary <- S7::new_class("CritStationary", parent = criterion,
  properties = list(tol = S7::class_numeric))

S7::method(crit_needs, CritStationary) <- function(criterion) "stationarity"

S7::method(crit_met, CritStationary) <- function(criterion, state) {
  s <- state$stationarity
  if (is.null(s) || !length(s) || !is.finite(s)) return(FALSE)
  s < criterion@tol
}

#' @title Stop When the Method's Own Measure of Progress Is Small
#'
#' @description
#' The stopping rule for a method that has no gradient to test.
#'
#' @param tol Numeric tolerance. Defaults to \code{1e-8}.
#'
#' @details
#' A gradient-based method detects its arrival through the vanishing of
#' \eqn{\nabla f}.
#' None of the derivative-free methods can use that test, and for the
#' non-smooth problems they exist to solve it would not be the right test even
#' if they could: at the minimum of \eqn{\lvert x \rvert} any evaluated
#' subgradient is \eqn{\pm 1}, so \code{\link{crit_grad}} would never fire at
#' the solution itself.
#'
#' Each such method therefore reports a non-negative scalar of its own that goes
#' to zero as it converges, and this rule tests that. What the scalar
#' \emph{is} differs, deliberately, because the natural measure differs:
#' \describe{
#'   \item{\code{\link{nelder_mead}}}{the diameter of the simplex, so the
#'     tolerance is on the parameter scale.}
#'   \item{\code{\link{compass}}}{the poll size \eqn{\Delta}. This is the
#'     rule with a theorem behind it: the limit points of a pattern search with
#'     \eqn{\Delta \to 0} are Clarke stationary.}
#'   \item{\code{\link{bundle}}}{the optimality estimate
#'     \eqn{\lVert p \rVert^2 + \alpha}, which vanishes exactly when zero lies
#'     in the convex hull of the collected subgradients with no linearization
#'     error. Note that this is \emph{not} the predicted decrease, which
#'     carries a factor of the trust parameter and can therefore be driven to
#'     zero by that parameter shrinking rather than by the point becoming
#'     stationary.}
#' }
#' The measure appears in the trace as the \code{stationarity} column, so a run
#' can be read afterwards without knowing which method produced it.
#'
#' @return A \code{\link{criterion}} object.
#'
#' @examples
#' crit_stationary()
#' crit_stationary(1e-10)
#'
#' @seealso \code{\link{nelder_mead}}, \code{\link{compass}},
#'   \code{\link{bundle}}, \code{\link{crit_grad}}
#' @export
crit_stationary <- function(tol = 1e-8) {
  check_tol(tol)
  CritStationary(label = paste0("stationarity < ", format(tol)), tol = tol)
}


# --- run the budget ---------------------------------------------------------

#' @title S7 Class for the Empty Criterion
#' @description The class \code{\link{crit_never}} instantiates.
#' @return An S7 object inheriting from \code{\link{criterion}}.
#' @seealso \code{\link{crit_never}}
#' @keywords internal
CritNever <- S7::new_class("CritNever", parent = criterion)

S7::method(crit_met, CritNever) <- function(criterion, state) FALSE

#' @title Never Stop Early
#'
#' @description
#' The rule that never fires, so a run ends only when it exhausts its iteration
#' budget.
#'
#' @details
#' This is not a placeholder. For a stochastic method there is often nothing
#' left to test: every quantity a convergence rule could look at — the
#' objective, the gradient — is a noisy estimate drawn from whichever
#' observations happened to be sampled, and a tolerance applied to one of those
#' measures the noise rather than the progress. Such a run is meant to be
#' governed by its budget, and saying so with an object is better than leaving a
#' real criterion in place that quietly never fires.
#'
#' A run that ends this way reports \code{converged = FALSE}, which is the
#' truth: the budget ran out, and nothing checked whether the answer was any
#' good. It is the same discipline everywhere else in the package — convergence
#' is what a rule confirmed, never what the run merely stopped doing.
#'
#' @return A \code{\link{criterion}} object.
#'
#' @examples
#' crit_never()
#'
#' @seealso \code{\link{adam}}, \code{\link{crit_grad}}
#' @export
crit_never <- function() CritNever(label = "iteration budget")


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

#' The Shared Body of the Two Combinators
#'
#' @description
#' Validates the arguments and builds the combined criterion, so that
#' \code{\link{crit_any}} and \code{\link{crit_all}} reject the same nonsense
#' in the same words.
#'
#' @param dots A list of \code{\link{criterion}} objects.
#' @param how Either \code{"any"} or \code{"all"}.
#'
#' @return A \code{\link{criterion}} object.
#'
#' @keywords internal
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
#' @details
#' \strong{The default rule of the gradient methods.} \code{\link{gd}},
#' \code{\link{cg}}, \code{\link{bb}}, \code{\link{newton}},
#' \code{\link{bfgs}} and \code{\link{lbfgs}} default to
#' \code{crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())}: the point is
#' stationary, or the objective has stopped moving, or the parameters have.
#' Since a disjunction can only get weaker as terms are added, a run that ends
#' under this rule would have ended under a gradient rule alone at best later
#' and never earlier.
#'
#' What that buys and what it costs was measured over the package's own
#' \code{\link{test_problems}}, six methods on eight problems. Against the
#' gradient rule alone it converges on 44 of the 48 runs rather than 41, and
#' costs 19370 objective evaluations rather than 22299. The three it gains are
#' \code{cg} and \code{bb} on the non-smooth \code{abs_sum} and \code{gd} on
#' Beale, and none is lost.
#'
#' The cost is that a run stops sooner, so the point it reports is further from
#' the solution. Measured, 10 of the 48 end at a gradient more than a hundred
#' times larger, and the worst are runs that were reaching absurd precision
#' anyway: \code{bfgs} on Rosenbrock ends at 8.1e-06 rather than 4.4e-10, with
#' the objective 3.2e-13 above its minimum rather than 1.2e-21. On one run the
#' difference is real rather than cosmetic -- \code{cg} on \code{abs_sum},
#' where the objective ends 4.5e-02 above the minimum rather than 1.9e-03, and
#' the flag reads \code{TRUE} where it used to read \code{FALSE}. That is a
#' smooth method on a non-smooth problem, where the objective stalls far from
#' the solution and a rule that reads a stall cannot tell the two apart. A
#' caller who needs stationarity asks for it: \code{criterion = crit_grad()}.
#'
#' \code{\link{crit_rel_obj}} was in that default until 0.6.0 and is not any
#' more, because it never fired: measured over the same 48 runs, the rule with
#' it and the rule without it agree on every count, every evaluation and every
#' reported point. It remains available and useful where an objective's scale
#' is not known in advance.
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
