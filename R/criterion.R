#' @title S7 Class for Convergence Criteria
#'
#' @import S7
#' @description
#' The abstract parent of every stopping rule. An optimizer carries one, and
#' a caller may replace it, combine several, or write a new kind: a class
#' that inherits from this one and has a method for [crit_met()] is a
#' stopping rule, and every algorithm in the package consults it.
#'
#' @details
#' A rule is an object that implements one generic, so a rule written by the
#' caller is treated in the same way as a shipped one. A new rule needs no
#' compiled code, whereas a new algorithm or a new kind of line search needs
#' a branch in the compiled loop.
#'
#' # The two generics
#'
#' [crit_met()] answers whether the run should stop, given the state of the
#' iteration just completed. [crit_needs()] names the `state` components that
#' the rule reads, so that an optimizer unable to supply one of them rejects
#' the rule when the run starts. The base method for `crit_needs()` returns
#' `character()`, so a rule that reads only the objective needs no method.
#'
#' # Combining
#'
#' [crit_any()] and [crit_all()] combine rules and are themselves criteria,
#' so combinations nest and their labels nest with them. [gd()], [cg()],
#' [bb()], [newton()], [bfgs()], [lbfgs()] and [chebyshev()] default to a
#' disjunction of three rules.
#'
#' @param label A short character label, reported as `criterion_met` when the
#'   rule fires and shown when the optimizer carrying it is printed.
#'
#' @return An S7 object of class `criterion`, carrying `label`. The class is
#'   abstract, so every value is an object of one of its subclasses or of a
#'   subclass written by the caller.
#'
#' @examples
#' # The class is abstract, so it cannot be instantiated directly...
#' try(criterion(label = "mine"))
#'
#' # ...but a class that inherits from it is a criterion, including a rule
#' # written by the caller.
#' Tiny <- S7::new_class("Tiny", parent = criterion,
#'                       properties = list(tol = S7::class_numeric))
#' S7::method(crit_met, Tiny) <- function(criterion, state)
#'   state$f_new < criterion@tol
#' crit_met(Tiny(label = "f < 1e-6", tol = 1e-6), list(f_new = 1e-9))
#'
#' @seealso [crit_grad()], [crit_rel_obj()],
#'   [crit_any()], [crit_met()]
#' @export
criterion <- S7::new_class(
  "criterion",
  properties = list(label = S7::class_character),
  abstract = TRUE
)


#' @title Has the Stopping Rule Been Met?
#'
#' @description
#' Evaluates a criterion on the state of the iteration just completed and
#' returns whether the run should stop. This is the one generic that a
#' criterion must implement: a class that inherits from [criterion()] and has
#' a method here is a stopping rule, and every algorithm in the package
#' consults it.
#'
#' @param criterion A [criterion()] object.
#' @param state A named list describing the current iteration; see Details.
#'
#' @details
#' `state` has the following components:
#' \describe{
#'   \item{`iter`}{the iteration just completed.}
#'   \item{`f_new`, `f_old`}{the objective after and before it.}
#'   \item{`x_new`, `x_old`}{the parameter vectors, likewise.}
#'   \item{`gradient`}{the gradient at `x_new`, or `NULL` when the
#'     method does not compute one. [prox_grad()] supplies the proximal
#'     gradient mapping instead.}
#'   \item{`stationarity`}{a non-negative measure of remaining progress,
#'     supplied in place of a gradient by [nelder_mead()], [compass()],
#'     [brent()], [sa()] and [bundle()], or `NULL`. See [crit_stationary()].}
#' }
#' The shipped methods supply `f_old` and `x_old` from the first call on,
#' except [brent()], which supplies them from its second iteration. A descent
#' method calls the rule with `f_old` and `x_old` set to `NULL` only after a
#' line search that accepted no step, so that only the state at `x_new` can
#' end the run there. With box constraints, `x_new`,
#' `x_old` and `gradient` are on the unconstrained scale on which the
#' optimizer moves.
#'
#' A rule that needs a component that a method does not supply, such as a
#' gradient on a derivative-free method, declares it through
#' [crit_needs()], and the optimizer then rejects the rule when the run
#' starts.
#'
#' @return A single logical.
#'
#' @examples
#' st <- list(iter = 3, f_new = 1.0000001, f_old = 1.0000002,
#'            x_new = c(1, 2), x_old = c(1, 2), gradient = c(1e-9, -2e-9))
#' crit_met(crit_grad(1e-8), st)
#' crit_met(crit_abs_obj(1e-12), st)
#'
#' @seealso [crit_needs()], [check_criterion()], [criterion()]
#' @export
crit_met <- S7::new_generic("crit_met", "criterion",
                            function(criterion, state) S7::S7_dispatch())


#' @title What a Criterion Needs From the Iteration
#'
#' @description
#' Returns the names of the `state` components that a criterion requires, so
#' that an algorithm can reject a rule that it cannot evaluate.
#'
#' @details
#' A derivative-free method computes no gradient, so [crit_grad()] given to
#' such a method would test `NULL` at every iteration and never stop the run.
#' [check_criterion()] compares the value of this generic with the value of
#' [optimizer_provides()] and, when the run starts, signals an error that
#' names the missing component and the method.
#'
#' The base method returns `character()`, so a rule that reads only the
#' objective needs no method here: every optimizer evaluates the objective,
#' and such a rule is accepted by every method. The shipped rules declare
#' two names, `"gradient"` and `"stationarity"`.
#'
#' @param criterion A [criterion()] object.
#'
#' @return A character vector of `state` component names, possibly empty.
#'
#' @examples
#' # A rule reading only the objective declares nothing; a gradient rule
#' # declares the one component it reads.
#' crit_needs(crit_abs_obj())
#' crit_needs(crit_grad())
#'
#' # A combination declares the union, so a method that lacks any one of
#' # them rejects it.
#' crit_needs(crit_any(crit_grad(), crit_stationary()))
#'
#' @seealso [crit_met()] for the rule itself, [check_criterion()] for the
#'   rejection this feeds, [optimizer_provides()] for the other half of the
#'   comparison.
#' @export
#' @aliases crit_needs.criterion
crit_needs <- S7::new_generic("crit_needs", "criterion",
                              function(criterion) S7::S7_dispatch())

S7::method(crit_needs, criterion) <- function(criterion) character()


# --- gradient ---------------------------------------------------------------

#' @title S7 Class for the Gradient Criterion
#'
#' @description
#' The class of the rule built by [crit_grad()], and its two methods. The
#' rule reads `state$gradient`, takes its max-norm or its 2-norm, and fires
#' when that norm falls below `tol`. It declares `"gradient"` through
#' [crit_needs()], so a derivative-free optimizer rejects it when the run
#' starts.
#'
#' @details
#' `crit_met()` returns `FALSE` when the gradient is `NULL`, empty or
#' contains an `NA`, so a state without a gradient never satisfies the rule.
#' The `norm` property selects the norm: on a gradient of \eqn{(3, 4)},
#' `crit_grad(4.5, "max")` fires and `crit_grad(4.5, "2")` does not.
#'
#' @param tol The tolerance, a single positive number.
#' @param norm Either `"max"` or `"2"`.
#'
#' @return An S7 object of class `CritGrad` inheriting from [criterion()],
#'   carrying `label`, `tol` and `norm`.
#'
#' @examples
#' st <- list(iter = 1, f_new = 1, f_old = 2, x_new = 1, x_old = 0,
#'            gradient = c(3, 4))
#' crit_needs(crit_grad())
#' c(max = crit_met(crit_grad(4.5, "max"), st),
#'   two = crit_met(crit_grad(4.5, "2"), st))
#'
#' # A gradient that is not there is not a small gradient.
#' crit_met(crit_grad(), list(f_new = 1, gradient = NULL))
#'
#' @seealso [crit_grad()] for the constructor and the tolerance it can
#'   attain, [crit_stationary()] for what a derivative-free method reads.
#' @name CritGrad-class
#' @aliases CritGrad crit_met.CritGrad crit_needs.CritGrad
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
#' Builds the rule \eqn{\lVert \nabla f \rVert < \texttt{tol}}, the test for
#' a stationary point and the rule for a maximum likelihood fit. Only a
#' method that computes a gradient can use it; a derivative-free optimizer
#' rejects it when the run starts.
#'
#' @param tol Numeric tolerance, a single positive number. Defaults to
#'   `1e-6`; see Details for the reason.
#' @param norm `"max"` (default) or `"2"`. Partial matching applies, and any
#'   other string is rejected by [match.arg()].
#'
#' @details
#' The max-norm is the default because it does not grow with the dimension,
#' as the 2-norm does: the same tolerance then has the same meaning for a
#' problem with two parameters and for one with two hundred, whereas `1e-6`
#' in the 2-norm is a much stricter condition in high dimension.
#'
#' The smallest gradient that a run can reach is set by the objective and not
#' by the method. A line search accepts a step only when the objective
#' decreases by a definite amount, and near a minimum that decrease is about
#' \eqn{\lVert \nabla f \rVert^{2} / (2\lambda)} for a curvature
#' \eqn{\lambda}. Once it falls below the rounding error of the objective,
#' about \eqn{\varepsilon \lvert f \rvert}, no step can be verified and the
#' search stops, so the smallest attainable gradient is about
#' \eqn{\sqrt{2 \lambda \varepsilon \lvert f^{*} \rvert}} and grows with the
#' value at the solution. Adding a constant to the objective moves neither
#' the minimizer nor the gradient, but it changes this bound (the last
#' example below). The default suits an objective of order one at its
#' solution, such as a log-likelihood divided by the number of observations;
#' an objective whose minimum is of the order of millions needs a larger
#' tolerance, and one whose minimum is zero admits a much smaller one.
#'
#' @return An S7 object of class [CritGrad], inheriting from [criterion()].
#'
#' @examples
#' crit_grad()
#' crit_grad(1e-10, norm = "2")
#'
#' # The two norms differ, and the max-norm is the looser of the two.
#' st <- list(f_new = 1, f_old = 2, x_new = 1, x_old = 0, gradient = c(3, 4))
#' c(max = crit_met(crit_grad(4.5, "max"), st),
#'   two = crit_met(crit_grad(4.5, "2"), st))
#'
#' # Asking for the rule alone, where the shipped default is a disjunction.
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#' r <- minimize(bfgs(criterion = crit_grad()), f, c(-1.2, 1), gr = g)
#' c(r@criterion_met, max(abs(g(r@par))))
#'
#' # And the constant added to the objective, which moves neither the
#' # minimizer nor the gradient, changes what is reachable.
#' shifted <- function(p) f(p) + 1e6
#' s <- minimize(cg(criterion = crit_grad(1e-14), maxit = 20000), shifted,
#'               c(-1.2, 1), gr = g)
#' c(converged = s@converged, gradient = max(abs(g(s@par))))
#'
#' @seealso [crit_stationary()] for the rule of a derivative-free method,
#'   [crit_any()] for the disjunction that the gradient methods use by
#'   default, [test_problems()] for the problems of the examples.
#' @export
crit_grad <- function(tol = 1e-6, norm = c("max", "2")) {
  norm <- match.arg(norm)
  check_tol(tol)
  CritGrad(label = paste0("gradient (", norm, "-norm) < ", format(tol)),
           tol = tol, norm = norm)
}


# --- objective --------------------------------------------------------------

#' @title S7 Class for the Absolute Objective Criterion
#'
#' @description
#' The class of the rule built by [crit_abs_obj()]. The rule reads
#' `state$f_new` and `state$f_old` and fires when they differ by less than
#' `tol`. It declares nothing through [crit_needs()], because every optimizer
#' evaluates the objective, so every method accepts it.
#'
#' @details
#' `crit_met()` returns `FALSE` when `f_old` is `NULL` or not finite, so the
#' rule cannot fire without two values to compare.
#'
#' The tolerance is in the units of the objective, so the rule depends on
#' the scale of the problem: `1e-10` has a different meaning for a
#' log-likelihood of order one and for one of order \eqn{10^{5}}.
#'
#' @param tol The tolerance, a single positive number.
#'
#' @return An S7 object of class `CritAbsObj` inheriting from [criterion()],
#'   carrying `label` and `tol`.
#'
#' @examples
#' crit_needs(crit_abs_obj())
#' crit_met(crit_abs_obj(1e-6),
#'          list(f_new = 1.0000001, f_old = 1.0000002))
#'
#' # Without a previous value there is nothing to compare.
#' crit_met(crit_abs_obj(), list(f_new = 1, f_old = NULL))
#'
#' @seealso [crit_abs_obj()] for the constructor, [crit_rel_obj()] for the
#'   scale-free version.
#' @name CritAbsObj-class
#' @aliases CritAbsObj crit_met.CritAbsObj
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
#' Builds the rule \eqn{\lvert f_{new} - f_{old} \rvert < \texttt{tol}}. It
#' reads only the objective, so every optimizer accepts it, and it is one of
#' the three terms of the default rule of the gradient methods.
#'
#' @details
#' The tolerance is in the units of the objective, so it depends on the
#' scale of the problem: `1e-10` is a strict condition on a log-likelihood
#' divided by the number of observations and a loose one on a summed
#' log-likelihood of order \eqn{10^{5}}. [crit_rel_obj()] is the relative
#' version.
#'
#' A rule on the objective cannot distinguish a stalled run from a converged
#' one. Where that distinction matters, use [crit_grad()].
#'
#' @param tol Numeric tolerance, a single positive number. Defaults to
#'   `1e-10`.
#'
#' @return An S7 object of class [CritAbsObj], inheriting from
#'   [criterion()].
#'
#' @examples
#' crit_abs_obj()
#' crit_abs_obj(1e-6)
#'
#' crit_met(crit_abs_obj(1e-6), list(f_new = 1.0000001, f_old = 1.0000002))
#'
#' @seealso [crit_rel_obj()] for the relative version, [crit_grad()] for
#'   the test of stationarity, [crit_any()] for the default rule that
#'   contains this one.
#' @export
crit_abs_obj <- function(tol = 1e-10) {
  check_tol(tol)
  CritAbsObj(label = paste0("|df| < ", format(tol)), tol = tol)
}


#' @title S7 Class for the Relative Objective Criterion
#'
#' @description
#' The class of the rule built by [crit_rel_obj()]. The rule fires when
#' \eqn{\lvert f_{new} - f_{old}\rvert < \texttt{tol}\,
#' (\lvert f_{old}\rvert + \texttt{tol})}, so the change is compared with the
#' size of the objective and one tolerance serves for any units. It declares
#' nothing through [crit_needs()].
#'
#' @details
#' `crit_met()` returns `FALSE` when `f_old` is `NULL` or not finite.
#'
#' The term `+ tol` is a floor on the scale. Without it the threshold would
#' be zero at an objective value of zero, and since the comparison is strict
#' the rule could never fire at an optimum of zero.
#'
#' @param tol The tolerance, a single positive number.
#'
#' @return An S7 object of class `CritRelObj` inheriting from [criterion()],
#'   carrying `label` and `tol`.
#'
#' @examples
#' crit_met(crit_rel_obj(1e-6), list(f_new = 1.0000001, f_old = 1.0000002))
#'
#' # With the floor, the rule can fire at an optimum of zero.
#' crit_met(crit_rel_obj(), list(f_new = 1e-30, f_old = 0))
#'
#' @seealso [crit_rel_obj()] for the constructor, [crit_abs_obj()] for the
#'   version in the units of the objective.
#' @name CritRelObj-class
#' @aliases CritRelObj crit_met.CritRelObj
#' @keywords internal
CritRelObj <- S7::new_class("CritRelObj", parent = criterion,
  properties = list(tol = S7::class_numeric))

S7::method(crit_met, CritRelObj) <- function(criterion, state) {
  if (is.null(state$f_old) || !is.finite(state$f_old)) return(FALSE)
  # The floor keeps the threshold positive at an optimum of zero, where the
  # strict comparison could otherwise never hold.
  abs(state$f_new - state$f_old) <
    criterion@tol * (abs(state$f_old) + criterion@tol)
}

#' @title Stop When the Objective Stops Moving (Relative)
#'
#' @description
#' The rule
#' \eqn{\lvert f_{new} - f_{old} \rvert < \texttt{tol}\,(\lvert f_{old} \rvert + \texttt{tol})}.
#'
#' @details
#' The term `+ tol` is a floor on the scale. Without it the threshold would
#' be zero at an objective value of zero, and since the comparison is strict
#' the rule could never fire at an optimum of zero.
#'
#' The rule was part of the default rule of the gradient methods until
#' version 0.6.0. It was removed because on the package's [test_problems()]
#' it never fired before the other terms. It is useful where the scale of the
#' objective is not known in advance, which is where [crit_abs_obj()] is
#' hard to set.
#'
#' @param tol Numeric tolerance, a single positive number. Defaults to
#'   `1e-12`, smaller than the default of [crit_abs_obj()] because the
#'   quantity compared is a ratio.
#'
#' @return An S7 object of class [CritRelObj], inheriting from
#'   [criterion()].
#'
#' @examples
#' crit_rel_obj()
#'
#' # The same absolute change, at two objective scales.
#' c(small = crit_met(crit_rel_obj(1e-6), list(f_new = 1, f_old = 1 + 1e-7)),
#'   large = crit_met(crit_rel_obj(1e-6),
#'                    list(f_new = 1e6, f_old = 1e6 + 1e-7)))
#'
#' # With the floor, the rule can fire at an optimum of zero.
#' crit_met(crit_rel_obj(), list(f_new = 1e-30, f_old = 0))
#'
#' @seealso [crit_abs_obj()] for the version in the units of the objective,
#'   [crit_any()] for the default rule of the gradient methods.
#' @export
crit_rel_obj <- function(tol = 1e-12) {
  check_tol(tol)
  CritRelObj(label = paste0("|df| < ", format(tol), " (relative)"), tol = tol)
}


# --- parameters -------------------------------------------------------------

#' @title S7 Class for the Absolute Parameter Criterion
#'
#' @description
#' The class of the rule built by [crit_abs_par()]. The rule fires when the
#' largest change of a coordinate,
#' \eqn{\max_j \lvert x_j^{new} - x_j^{old}\rvert}, falls below `tol`, so the
#' tolerance is in the units of the parameters. It declares nothing through
#' [crit_needs()].
#'
#' @details
#' `crit_met()` returns `FALSE` when `x_old` is `NULL`. With box constraints
#' the comparison is on the **unconstrained** scale, on which the optimizer
#' moves, so the same tolerance has different meanings for a variance near
#' zero and for one near a thousand.
#'
#' @param tol The tolerance, a single positive number.
#'
#' @return An S7 object of class `CritAbsPar` inheriting from [criterion()],
#'   carrying `label` and `tol`.
#'
#' @examples
#' crit_met(crit_abs_par(1e-6),
#'          list(x_new = c(1, 2), x_old = c(1, 2 + 1e-9)))
#' crit_met(crit_abs_par(), list(x_new = c(1, 2), x_old = NULL))
#'
#' @seealso [crit_abs_par()] for the constructor, [crit_rel_par()] for the
#'   version relative to each coordinate.
#' @name CritAbsPar-class
#' @aliases CritAbsPar crit_met.CritAbsPar
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
#' Builds the rule
#' \eqn{\max_j \lvert x_j^{new} - x_j^{old} \rvert < \texttt{tol}}, so the
#' tolerance is in the units of the parameters. It reads only the iterates,
#' and it is one of the three terms of the default rule of the gradient
#' methods.
#'
#' @details
#' With box constraints the comparison is on the **unconstrained** scale, on
#' which the optimizer moves. Near a bound a short distance in \eqn{\theta}
#' corresponds to a long distance in \eqn{\eta}, so the rule fires later
#' there than the change on the parameter scale suggests.
#'
#' A method whose iterate does not move at an iteration, such as [brent()]
#' after a trial point that does not improve on the best one, reports a
#' change of zero, so the rule fires at that iteration.
#'
#' @param tol Numeric tolerance, a single positive number. Defaults to
#'   `1e-8`.
#'
#' @return An S7 object of class [CritAbsPar], inheriting from
#'   [criterion()].
#'
#' @examples
#' crit_abs_par()
#' crit_met(crit_abs_par(1e-6), list(x_new = c(1, 2), x_old = c(1, 2 + 1e-9)))
#'
#' @seealso [crit_rel_par()] for the version relative to each coordinate,
#'   [crit_any()] for the default rule that contains this one.
#' @export
crit_abs_par <- function(tol = 1e-8) {
  check_tol(tol)
  CritAbsPar(label = paste0("|dx| < ", format(tol)), tol = tol)
}


#' @title S7 Class for the Relative Parameter Criterion
#'
#' @description
#' The class of the rule built by [crit_rel_par()]. The rule fires when the
#' change of every coordinate, divided by the size of that coordinate, falls
#' below `tol`, so a parameter of order \eqn{10^{3}} and one of order
#' \eqn{10^{-3}} are held to the same number of digits. It declares nothing
#' through [crit_needs()].
#'
#' @details
#' The test is
#' \eqn{\max_j \lvert x_j^{new} - x_j^{old}\rvert /
#' (\lvert x_j^{old}\rvert + \texttt{tol}) < \texttt{tol}}, with the floor of
#' [CritRelObj], which keeps the denominator positive at a coordinate of
#' zero. `crit_met()` returns `FALSE` when `x_old` is `NULL`.
#'
#' @param tol The tolerance, a single positive number, used both as the floor
#'   and as the threshold.
#'
#' @return An S7 object of class `CritRelPar` inheriting from [criterion()],
#'   carrying `label` and `tol`.
#'
#' @examples
#' # The same relative change at two very different scales.
#' crit_met(crit_rel_par(1e-6), list(x_new = 1000, x_old = 1000.0001))
#' crit_met(crit_rel_par(1e-6), list(x_new = 1e-3, x_old = 1e-3 + 1e-10))
#'
#' @seealso [crit_rel_par()] for the constructor, [crit_abs_par()] for the
#'   version in the units of the parameters.
#' @name CritRelPar-class
#' @aliases CritRelPar crit_met.CritRelPar
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
#' @details
#' The term `+ tol` in each denominator is a floor, as in [crit_rel_obj()],
#' which keeps the denominator positive at a coordinate of zero. The same
#' number therefore serves as the floor and as the threshold.
#'
#' @param tol Numeric tolerance, a single positive number. Defaults to
#'   `1e-8`.
#'
#' @return An S7 object of class [CritRelPar], inheriting from
#'   [criterion()].
#'
#' @examples
#' crit_rel_par()
#'
#' # The same relative change at two very different scales.
#' c(large = crit_met(crit_rel_par(1e-6), list(x_new = 1000,
#'                                             x_old = 1000.0001)),
#'   small = crit_met(crit_rel_par(1e-6), list(x_new = 1e-3,
#'                                             x_old = 1e-3 + 1e-10)))
#'
#' @seealso [crit_abs_par()] for the version in the units of the parameters.
#' @export
crit_rel_par <- function(tol = 1e-8) {
  check_tol(tol)
  CritRelPar(label = paste0("|dx| < ", format(tol), " (relative)"), tol = tol)
}


# --- stationarity -----------------------------------------------------------

#' @title S7 Class for the Stationarity Criterion
#'
#' @description
#' The class of the rule built by [crit_stationary()], and its two methods.
#' The rule reads `state$stationarity`, the non-negative measure that a
#' method without a gradient test reports instead, and fires when it falls
#' below `tol`. It declares `"stationarity"` through [crit_needs()], so a
#' method that reports no such measure rejects it when the run starts.
#'
#' @details
#' `crit_met()` returns `FALSE` when the measure is `NULL`, empty or not
#' finite. The measure differs by method, so the tolerance has a different
#' meaning for each: the largest distance from the best vertex to another
#' vertex of the simplex for [nelder_mead()], the poll size for [compass()],
#' the half-width of the bracket for [brent()], Corana's termination measure
#' for [sa()], and the optimality estimate \eqn{\lVert p\rVert^{2} + \alpha}
#' for [bundle()].
#'
#' @param tol The tolerance, a single positive number.
#'
#' @return An S7 object of class `CritStationary` inheriting from
#'   [criterion()], carrying `label` and `tol`.
#'
#' @examples
#' crit_needs(crit_stationary())
#' crit_met(crit_stationary(1e-6), list(stationarity = 1e-9))
#' crit_met(crit_stationary(), list(stationarity = NULL))
#'
#' @seealso [crit_stationary()] for the constructor and the measure of each
#'   method, [crit_grad()] for the gradient-based rule.
#' @name CritStationary-class
#' @aliases CritStationary crit_met.CritStationary crit_needs.CritStationary
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
#' Builds the stopping rule for a method that has no gradient to test. Each
#' such method reports a non-negative scalar of its own that tends to zero as
#' the method converges, and this rule tests that scalar. It is the default
#' of [nelder_mead()], [compass()], [brent()], [sa()] and [bundle()]. The
#' gradient-based methods report no such scalar, so they reject the rule
#' when the run starts.
#'
#' @param tol Numeric tolerance, a single positive number. Defaults to
#'   `1e-8`. The quantity to which the tolerance applies depends on the
#'   method; see Details.
#'
#' @details
#' A gradient-based method detects convergence through the vanishing of
#' \eqn{\nabla f}. The derivative-free methods cannot use that test, and on
#' a non-smooth problem it is not a reliable test in any case: at every point
#' other than the minimum of \eqn{\lvert x \rvert} the gradient is
#' \eqn{\pm 1}, so [crit_grad()] fires only at an iterate that lands exactly
#' on the minimum.
#'
#' The scalar differs by method:
#' \describe{
#'   \item{[nelder_mead()]}{the largest distance from the best vertex to
#'     another vertex of the simplex, so the tolerance is on the parameter
#'     scale.}
#'   \item{[compass()]}{the poll size \eqn{\Delta}. For a continuously
#'     differentiable objective, the limit points of a coordinate pattern
#'     search with \eqn{\Delta \to 0} are stationary (Torczon 1997).}
#'   \item{[brent()]}{the half-width of the bracket.}
#'   \item{[sa()]}{Corana's termination measure; see [sa()].}
#'   \item{[bundle()]}{the optimality estimate
#'     \eqn{\lVert p \rVert^2 + \alpha}, which vanishes exactly when zero lies
#'     in the convex hull of the collected subgradients with no linearization
#'     error. It is not the predicted decrease, which contains a factor of the
#'     trust parameter and therefore tends to zero when that parameter
#'     shrinks, whether or not the point is stationary.}
#' }
#' The measure appears in the trace as the `stationarity` column.
#'
#' @return An S7 object of class [CritStationary], inheriting from
#'   [criterion()].
#'
#' @examples
#' crit_stationary()
#' crit_stationary(1e-10)
#'
#' # It reads state$stationarity, and an absent measure is not a small one.
#' c(present = crit_met(crit_stationary(1e-6), list(stationarity = 1e-9)),
#'   absent  = crit_met(crit_stationary(), list(stationarity = NULL)))
#'
#' # The measure differs by method (a size of the simplex for nelder_mead(),
#' # the poll size for compass()), and the trace records it in the same
#' # column.
#' q <- function(p) sum((p - c(1, 2))^2)
#' names(minimize(nelder_mead(keep_trace = TRUE), q, c(0, 0))@trace)
#' names(minimize(compass(keep_trace = TRUE), q, c(0, 0))@trace)
#'
#' @seealso [nelder_mead()], [compass()], [brent()], [sa()] and [bundle()]
#'   for the five measures, [crit_grad()] for the gradient-based rule.
#' @export
crit_stationary <- function(tol = 1e-8) {
  check_tol(tol)
  CritStationary(label = paste0("stationarity < ", format(tol)), tol = tol)
}


# --- run the budget ---------------------------------------------------------

#' @title S7 Class for the Empty Criterion
#'
#' @description
#' The class of the rule built by [crit_never()]. Its `crit_met()` returns
#' `FALSE` at every state without reading it, so a run with this rule ends
#' when a budget runs out or when the method stops for a reason of its own.
#' It declares nothing through [crit_needs()], so every optimizer accepts it.
#'
#' @details
#' It has no `tol`, since it compares nothing, and it is the only criterion
#' class with no property beyond `label`. [adam()] is the only shipped
#' method that uses it by default.
#'
#' @return An S7 object of class `CritNever` inheriting from [criterion()],
#'   carrying `label` alone.
#'
#' @examples
#' crit_needs(crit_never())
#' crit_met(crit_never(), list(f_new = 0, f_old = 0, gradient = c(0, 0)))
#'
#' @seealso [crit_never()] for the constructor, [adam()] for the method that
#'   uses it.
#' @name CritNever-class
#' @aliases CritNever crit_met.CritNever
#' @keywords internal
CritNever <- S7::new_class("CritNever", parent = criterion)

S7::method(crit_met, CritNever) <- function(criterion, state) FALSE

#' @title The Rule That Never Fires
#'
#' @description
#' Builds the rule that never fires, so a run ends when it exhausts a budget
#' and reports `converged = FALSE`. It is the default of [adam()], a method
#' for stochastic objectives, on which every quantity that a stopping rule
#' could read is an estimate.
#'
#' @details
#' On a stochastic objective the value and the gradient are noisy estimates
#' computed from the observations that were sampled, and a tolerance applied
#' to one of them tests the noise. Such a run is governed by its budget, and
#' this rule states that explicitly.
#'
#' A run that ends on its budget reports `converged = FALSE`, because no rule
#' tested the answer. A descent method can also stop before its budget, when
#' its line search accepts no step, and it then reports `converged = FALSE`
#' as well. In the package `converged` is `TRUE` only when a stopping rule
#' fired, when a line search with a declared resolution established that no
#' step decreases the objective by more than that resolution, or when
#' [brent()] reduced its bracket to the resolution of the iterate.
#'
#' @return An S7 object of class [CritNever], inheriting from [criterion()].
#'
#' @examples
#' crit_never()
#'
#' # It fires at no state, even at an exact minimum.
#' crit_met(crit_never(), list(f_new = 0, f_old = 0, x_new = 1, x_old = 1,
#'                             gradient = c(0, 0)))
#'
#' # A run with this rule ends on its budget, and the result reports it.
#' r <- minimize(adam(maxit = 50), function(p) sum(p^2), c(1, 1),
#'               gr = function(p) 2 * p)
#' c(r@converged, r@criterion_met)
#'
#' @seealso [adam()] for the method that defaults to it, [crit_grad()] for
#'   the rule to pass instead when the objective is exact.
#' @export
crit_never <- function() CritNever(label = "iteration budget")


# --- combinators ------------------------------------------------------------

#' @title S7 Class for a Combination of Criteria
#'
#' @description
#' The class built by both [crit_any()] and [crit_all()], and its two
#' methods. Its `crit_met()` evaluates every rule that it contains and
#' reduces the results with `any()` or `all()` according to `how`; its
#' `crit_needs()` is the union of the needs of those rules, so a derivative-free
#' method rejects a combination that contains a gradient rule, as it rejects
#' the gradient rule alone.
#'
#' @details
#' A `CritCombine` is itself a [criterion()], so combinations nest and the
#' label nests with them: `crit_any(crit_all(a, b), c)` reads
#' `a and b or c`. Every rule is evaluated at every call, because `any()` and
#' `all()` receive the whole vector of results and do not short-circuit.
#'
#' @param criteria A list of [criterion()] objects.
#' @param how Either `"any"` or `"all"`.
#'
#' @return An S7 object of class `CritCombine` inheriting from [criterion()],
#'   carrying `label`, `criteria` and `how`.
#'
#' @examples
#' # The needs are the union, so a simplex method rejects this combination.
#' crit_needs(crit_any(crit_grad(), crit_stationary()))
#'
#' st <- list(iter = 3, f_new = 1, f_old = 2, x_new = 1, x_old = 1,
#'            gradient = c(1e-9, -2e-9))
#' c(any = crit_met(crit_any(crit_grad(1e-8), crit_never()), st),
#'   all = crit_met(crit_all(crit_grad(1e-8), crit_never()), st))
#'
#' # Combinations nest, and so does the label.
#' crit_any(crit_all(crit_grad(), crit_abs_par()), crit_never())@label
#'
#' @seealso [crit_any()] and [crit_all()] for the constructors,
#'   [combine_criteria()] for the shared body.
#' @name CritCombine-class
#' @aliases CritCombine crit_met.CritCombine crit_needs.CritCombine
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
#' [crit_any()] and [crit_all()] reject the same invalid input with the same
#' messages: an empty call and an argument that is not a criterion.
#'
#' @details
#' The label is the sub-labels joined by ` or ` or ` and ` according to
#' `how`, which is why a nested combination reads as one sentence.
#'
#' @param dots A list of [criterion()] objects.
#' @param how Either `"any"` or `"all"`.
#'
#' @return An S7 object of class [CritCombine], inheriting from
#'   [criterion()].
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
#' Combines criteria disjunctively: the run ends as soon as any of the rules
#' is satisfied.
#'
#' @details
#' **The default rule of the gradient methods.** [gd()], [cg()], [bb()],
#' [newton()], [bfgs()], [lbfgs()] and [chebyshev()] default to
#' `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`: the point is
#' stationary, or the objective has stopped moving, or the parameters have.
#' A disjunction can only become weaker as terms are added, so a run that
#' ends under this rule would have ended under the gradient rule alone at
#' the same iteration or later, never earlier.
#'
#' On the package's [test_problems()], with the six methods other than
#' [chebyshev()], this rule converges on every run on which the gradient rule
#' alone converges, and on three more (`cg` and `bb` on `abs_sum`, `gd` on
#' `beale`), and it needs fewer evaluations of the objective and the gradient
#' in total. The cost is that a run stops sooner, at a point further from the
#' solution, usually by a margin far below the statistical precision of an
#' estimate. On a non-smooth problem a smooth method can stall far from the
#' solution, and a rule that reads a stall reports it as convergence (`cg` on
#' `abs_sum` is such a case). A caller who needs a stationary point uses
#' `criterion = crit_grad()`.
#'
#' [crit_rel_obj()] was part of this default until 0.6.0. It was removed
#' because on the same problems it never fired before the other terms.
#'
#' @param ... [criterion()] objects.
#'
#' @return A [criterion()] object, so combinations nest.
#'
#' @examples
#' crit_any(crit_grad(1e-8), crit_rel_obj(1e-12))
#'
#' @seealso [crit_all()]
#' @export
crit_any <- function(...) combine_criteria(list(...), "any")

#' @title Stop Only When Every Rule Fires
#'
#' @description
#' Combines criteria conjunctively, for a run that should stop only when
#' several rules hold at once. A conjunction can only become **stronger** as
#' terms are added, so a run that ends under it would have ended under any
#' one of its terms alone at the same iteration or earlier.
#'
#' @details
#' It is useful when one rule alone is known to fire early:
#' `crit_all(crit_grad(), crit_abs_par())` requires a stationary point at
#' which the iterate has also stopped moving, a condition that a run circling
#' a flat optimum does not satisfy.
#'
#' The needs of a conjunction are the union of the needs of its rules, so an
#' optimizer rejects a conjunction that contains a rule that it cannot
#' evaluate, as it rejects that rule alone.
#'
#' @param ... [criterion()] objects. At least one is required; anything that
#'   is not a criterion raises an error.
#'
#' @return An S7 object of class [CritCombine], inheriting from
#'   [criterion()], so combinations nest.
#'
#' @examples
#' crit_all(crit_grad(1e-6), crit_abs_par(1e-10))
#'
#' # Both must hold. Here the gradient rule fires and the parameter rule does
#' # not, so the conjunction does not.
#' st <- list(f_new = 1, f_old = 2, x_new = c(1, 2), x_old = c(1, 3),
#'            gradient = c(1e-9, -2e-9))
#' c(grad = crit_met(crit_grad(1e-8), st),
#'   par  = crit_met(crit_abs_par(1e-10), st),
#'   all  = crit_met(crit_all(crit_grad(1e-8), crit_abs_par(1e-10)), st))
#'
#' @seealso [crit_any()] for the disjunction, [criterion()] for what a rule
#'   is.
#' @export
crit_all <- function(...) combine_criteria(list(...), "all")


#' @title Print Method for Criteria
#' @name print.criterion
#'
#' @description
#' Shows a criterion in one line, as `<criterion>` followed by its label.
#' When the rule ends a run, the result reports the same label in
#' `criterion_met`.
#'
#' @param x A [criterion()] object.
#' @param ... Unused.
#'
#' @return `x`, invisibly. Called for the output.
#'
#' @examples
#' crit_grad()
#' crit_stationary(1e-10)
#'
#' # A combination prints as one sentence, and nests.
#' crit_any(crit_grad(), crit_rel_obj())
#' crit_any(crit_all(crit_grad(), crit_abs_par()), crit_never())
#'
#' # The result of a run that ends on the rule reports the same string.
#' r <- minimize(bfgs(criterion = crit_grad()), function(p) sum(p^2), c(1, 1),
#'               gr = function(p) 2 * p)
#' r@criterion_met
#'
#' @keywords internal
S7::method(print, criterion) <- function(x, ...) {
  cat("<criterion> ", x@label, "\n", sep = "")
  invisible(x)
}


#' Validate a Tolerance
#'
#' @description
#' Checks that the value is a single positive number, so that every criterion
#' constructor rejects invalid values with the same message. Zero is rejected
#' as well as a negative value, because a rule with a tolerance of zero can
#' never fire.
#'
#' @details
#' It is also called by [bundle()] for `t0`, `t_min`, `t_max` and `qp_tol`,
#' by [bb()] for `alpha0`, `alpha_min`, `alpha_max` and `curv_tol`, by
#' [nelder_mead()] for `degenerate_tol` and by [multistart()] for `spread`
#' and `distinct_tol`. These callers pass the name of their own argument in
#' `nm`, so that the message names the argument that was supplied.
#'
#' @param tol The value supplied.
#' @param nm Its name, for the message. Defaults to `"tol"`.
#'
#' @return Invisibly `TRUE`. Raises an error naming `nm` otherwise.
#'
#' @keywords internal
check_tol <- function(tol, nm = "tol") {
  if (length(tol) != 1L || !is.numeric(tol) || is.na(tol) || tol <= 0) {
    stop("'", nm, "' must be a single positive number.", call. = FALSE)
  }
  invisible(TRUE)
}
