#' @title S7 Class for Convergence Criteria
#'
#' @import S7
#' @description
#' The abstract parent of every stopping rule. An optimizer carries one, and
#' a caller may replace it, combine several, or write a new kind: a class
#' inheriting from this one with a method for [crit_met()] is a stopping rule
#' and every algorithm in the package will consult it.
#'
#' @details
#' The alternative would be an argument taking a string and a `switch` inside
#' every algorithm, which fixes the set of rules at the moment the package is
#' written and lets nothing outside add to it. Here a rule is an object
#' implementing one generic, so a rule of your own is treated exactly as a
#' shipped one. This is the package's most open extension point: the
#' algorithms need a branch in compiled code and the line searches do too,
#' but a criterion needs neither.
#'
#' # The two generics
#'
#' [crit_met()] answers whether the run should stop, given the state of the
#' iteration just completed. [crit_needs()] names the `state` components the
#' rule reads, so that an optimizer unable to fill one in can refuse the rule
#' when the run starts instead of accepting one that never fires. The base
#' method for `crit_needs()` returns `character()`, so a rule reading only
#' the objective needs no method at all.
#'
#' # Combining
#'
#' [crit_any()] and [crit_all()] combine rules and are themselves criteria,
#' so combinations nest and their labels nest with them. The gradient methods
#' default to a disjunction of three.
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
#' # ...but anything inheriting from it is a criterion, including a rule the
#' # package never anticipated.
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
#' Asks a criterion whether the run should stop, given the state of the
#' iteration just completed. This is the one generic a criterion must
#' implement: a class inheriting from [criterion()] with a method here is a
#' stopping rule, and every algorithm in the package will consult it.
#'
#' @param criterion A [criterion()] object.
#' @param state A named list describing the current iteration; see Details.
#'
#' @details
#' `state` carries everything any rule could need:
#' \describe{
#'   \item{`iter`}{the iteration just completed.}
#'   \item{`f_new`, `f_old`}{the objective after and before it.}
#'   \item{`x_new`, `x_old`}{the parameter vectors, likewise.}
#'   \item{`gradient`}{the gradient at `x_new`, or `NULL` when the
#'     method does not compute one.}
#'   \item{`stationarity`}{a non-negative measure of remaining progress,
#'     supplied by the derivative-free methods in place of a gradient, or
#'     `NULL`. See [crit_stationary()].}
#' }
#' A rule that needs something absent from `state` — a gradient, from a
#' derivative-free method — must say so through [crit_needs()] rather
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
#' @seealso [as_objective()], [crit_needs()], [check_criterion()]
#' @export
crit_met <- S7::new_generic("crit_met", "criterion",
                            function(criterion, state) S7::S7_dispatch())


#' @title What a Criterion Needs From the Iteration
#'
#' @description
#' The names of the `state` components a criterion requires, so that an
#' algorithm can reject a rule it cannot evaluate instead of accepting one that
#' never fires.
#'
#' @details
#' A derivative-free method has no gradient, so [crit_grad()] handed to one
#' would sit there testing `NULL` at every iteration and quietly never stop
#' the run. [check_criterion()] compares what this reports against what
#' [optimizer_provides()] reports and refuses the pair when the run starts,
#' naming both.
#'
#' The base method returns `character()`, so a rule reading only the
#' objective needs no method here: every optimizer evaluates the objective
#' and no rule reading it is ever refused. The two names the shipped rules
#' declare are `"gradient"` and `"stationarity"`.
#'
#' @param criterion A [criterion()] object.
#'
#' @return A character vector of `state` component names, possibly empty.
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
#' The rule [crit_grad()] builds, and the two methods it implements. It
#' reads `state$gradient`, takes its max-norm or its 2-norm, and fires when
#' that falls below `tol`. It declares `"gradient"` through [crit_needs()],
#' so a derivative-free optimizer refuses it when the run starts.
#'
#' @details
#' `crit_met()` returns `FALSE` when the gradient is `NULL`, empty or carries
#' an `NA`, so a state that has not yet filled it in never accidentally
#' satisfies the rule. Which norm is used is the `norm` property: on a
#' gradient of \eqn{(3, 4)}, `crit_grad(4.5, "max")` fires and
#' `crit_grad(4.5, "2")` does not.
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
#' a stationary point and the one a maximum likelihood fit should ask for.
#' Usable only by a method that computes a gradient; a derivative-free
#' optimizer refuses it when the run starts instead of accepting a rule that
#' can never fire.
#'
#' @param tol Numeric tolerance, a single positive number. Defaults to
#'   `1e-6`; see below for why that rather than something smaller.
#' @param norm `"max"` (default) or `"2"`. Partial matching applies and any
#'   other string is refused, naming both.
#'
#' @details
#' The max-norm is the default because it does not grow with the dimension the
#' way the 2-norm does: the same tolerance then means the same thing for a
#' two-parameter problem and a two-hundred-parameter one, whereas `1e-6` in
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
#' gradient --- takes the attainable gradient from `1.9e-9` at
#' \eqn{f^{*} = 0} to `4.4e-8` at \eqn{f^{*} = 1} and `6.5e-5` at
#' \eqn{f^{*} = 10^{6}}. The default suits an objective of order one at its
#' solution, which is what a log-likelihood per observation is; an objective
#' that lands in the millions needs a correspondingly looser tolerance, and one
#' that lands at zero can be asked for much more.
#'
#' Only usable by a method that computes a gradient; a derivative-free optimizer
#' rejects it rather than accepting a rule that can never fire.
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
#' @seealso [crit_stationary()] for what a derivative-free method reads,
#'   [crit_any()] for the disjunction the gradient methods default to,
#'   [test_problems()] for the battery the attainable-gradient figures come
#'   from.
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
#' The rule [crit_abs_obj()] builds. It reads `state$f_new` and
#' `state$f_old` and fires when they differ by less than `tol`. It declares
#' nothing through [crit_needs()], every optimizer evaluating the objective,
#' so no method refuses it.
#'
#' @details
#' `crit_met()` returns `FALSE` when `f_old` is `NULL` or not finite, which
#' is the state at the first iteration, so the rule cannot fire before there
#' are two values to compare.
#'
#' The tolerance is in the objective's own units, so this rule carries the
#' scale of the problem with it: `1e-10` means something different for a
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
#' # Nothing to compare against at the first iteration.
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
#' asks nothing of the method beyond an objective, so every optimizer accepts
#' it, and it is one of the three terms in the gradient methods' default
#' rule.
#'
#' @details
#' The tolerance is in the objective's own units, so it carries the scale of
#' the problem: `1e-10` is a strict demand on a log-likelihood per
#' observation and a loose one on a summed log-likelihood of order
#' \eqn{10^{5}}. [crit_rel_obj()] is the version that does not.
#'
#' A rule reading the objective cannot tell a stalled run from a converged
#' one. Where that distinction matters, ask for [crit_grad()].
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
#' @seealso [crit_rel_obj()] for the scale-free version, [crit_grad()] for
#'   the test of stationarity, [crit_any()] for the default rule that holds
#'   this one.
#' @export
crit_abs_obj <- function(tol = 1e-10) {
  check_tol(tol)
  CritAbsObj(label = paste0("|df| < ", format(tol)), tol = tol)
}


#' @title S7 Class for the Relative Objective Criterion
#'
#' @description
#' The rule [crit_rel_obj()] builds. It fires when
#' \eqn{\lvert f_{new} - f_{old}\rvert < \texttt{tol}\,
#' (\lvert f_{old}\rvert + \texttt{tol})}, so the comparison is against the
#' objective's own scale and one tolerance serves whatever units the problem
#' is in. It declares nothing through [crit_needs()].
#'
#' @details
#' `crit_met()` returns `FALSE` when `f_old` is `NULL` or not finite, as at
#' the first iteration.
#'
#' The `+ tol` in the denominator is a floor and is load-bearing: an
#' objective whose optimum sits at zero would otherwise be compared against a
#' vanishing scale, and the rule would either never fire or fire at once.
#'
#' @param tol The tolerance, a single positive number.
#'
#' @return An S7 object of class `CritRelObj` inheriting from [criterion()],
#'   carrying `label` and `tol`.
#'
#' @examples
#' crit_met(crit_rel_obj(1e-6), list(f_new = 1.0000001, f_old = 1.0000002))
#'
#' # The floor is what keeps an optimum at zero usable.
#' crit_met(crit_rel_obj(), list(f_new = 1e-30, f_old = 0))
#'
#' @seealso [crit_rel_obj()] for the constructor, [crit_abs_obj()] for the
#'   version in the objective's own units.
#' @name CritRelObj-class
#' @aliases CritRelObj crit_met.CritRelObj
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
#' @param tol Numeric tolerance. Defaults to `1e-12`.
#'
#' @details
#' The `+ tol` in the denominator is a floor, and it is load-bearing: an
#' objective whose optimum sits at zero would otherwise be compared against a
#' vanishing scale, and the rule would either never fire or fire at once.
#'
#' It was in the gradient methods' default rule until version 0.6.0 and is
#' not any more, because it never fired there: measured over the package's
#' own [test_problems()], six methods on eight problems, the default with it
#' and the default without it agree on every flag, every evaluation count and
#' every reported point. It remains useful where an objective's scale is not
#' known in advance, which is exactly where [crit_abs_obj()] is hard to set.
#'
#' @param tol Numeric tolerance, a single positive number. Defaults to
#'   `1e-12`, tighter than [crit_abs_obj()]'s because the quantity compared
#'   is a ratio.
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
#' # And the floor, which keeps an optimum at zero usable.
#' crit_met(crit_rel_obj(), list(f_new = 1e-30, f_old = 0))
#'
#' @seealso [crit_abs_obj()] for the version in the objective's own units,
#'   [crit_any()] for the default rule it left.
#' @export
crit_rel_obj <- function(tol = 1e-12) {
  check_tol(tol)
  CritRelObj(label = paste0("|df| < ", format(tol), " (relative)"), tol = tol)
}


# --- parameters -------------------------------------------------------------

#' @title S7 Class for the Absolute Parameter Criterion
#'
#' @description
#' The rule [crit_abs_par()] builds. It fires when the largest coordinate
#' change \eqn{\max_j \lvert x_j^{new} - x_j^{old}\rvert} falls below `tol`,
#' so the tolerance is in the parameters' own units. It declares nothing
#' through [crit_needs()].
#'
#' @details
#' `crit_met()` returns `FALSE` when `x_old` is `NULL`, as at the first
#' iteration. With box constraints the comparison is on the **unconstrained**
#' scale, that being where the optimizer moves, so the same tolerance means
#' different things about a variance near zero and one near a thousand.
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
#'   version scaled by each coordinate.
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
#' tolerance is in the parameters' own units. It asks nothing of the method,
#' and it is one of the three terms in the gradient methods' default rule.
#'
#' @details
#' With box constraints the comparison is on the **unconstrained** scale,
#' that being where the optimizer moves. A coordinate approaching a bound
#' travels a long way in \eqn{\eta} for a short way in \eqn{\theta}, so this
#' rule fires later there than a reader of the parameter scale would expect.
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
#' @seealso [crit_rel_par()] for the version scaled by each coordinate,
#'   [crit_any()] for the default rule that holds this one.
#' @export
crit_abs_par <- function(tol = 1e-8) {
  check_tol(tol)
  CritAbsPar(label = paste0("|dx| < ", format(tol)), tol = tol)
}


#' @title S7 Class for the Relative Parameter Criterion
#'
#' @description
#' The rule [crit_rel_par()] builds. It fires when every coordinate's change,
#' divided by that coordinate's own size, falls below `tol`, so a parameter
#' of order \eqn{10^{3}} and one of order \eqn{10^{-3}} are held to the same
#' number of digits. It declares nothing through [crit_needs()].
#'
#' @details
#' The test is
#' \eqn{\max_j \lvert x_j^{new} - x_j^{old}\rvert /
#' (\lvert x_j^{old}\rvert + \texttt{tol}) < \texttt{tol}}, with the same
#' floor [CritRelObj] uses and for the same reason: a coordinate sitting at
#' zero would otherwise be divided by nothing. `crit_met()` returns `FALSE`
#' when `x_old` is `NULL`.
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
#'   version in the parameters' own units.
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
#' The `+ tol` in each denominator is a floor, as in [crit_rel_obj()]: a
#' coordinate sitting at zero would otherwise be divided by nothing. The same
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
#' @seealso [crit_abs_par()] for the version in the parameters' own units.
#' @export
crit_rel_par <- function(tol = 1e-8) {
  check_tol(tol)
  CritRelPar(label = paste0("|dx| < ", format(tol), " (relative)"), tol = tol)
}


# --- stationarity -----------------------------------------------------------

#' @title S7 Class for the Stationarity Criterion
#'
#' @description
#' The rule [crit_stationary()] builds, and the two methods it implements. It
#' reads `state$stationarity`, the non-negative measure a derivative-free
#' method reports in place of a gradient, and fires when it falls below
#' `tol`. It declares `"stationarity"` through [crit_needs()], so a method
#' that reports none refuses it when the run starts.
#'
#' @details
#' `crit_met()` returns `FALSE` when the measure is `NULL`, empty or not
#' finite. What the measure *is* differs by method, so the tolerance means
#' something different for each: the simplex diameter for [nelder_mead()],
#' the poll size for [compass()], Corana's termination measure for [sa()],
#' and the optimality estimate \eqn{\lVert p\rVert^{2} + \alpha} for
#' [bundle()]. [crit_stationary()] carries the comparison across all four.
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
#' @seealso [crit_stationary()] for the constructor and what each method
#'   reports, [crit_grad()] for the gradient-based rule.
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
#' Builds the stopping rule for a method that has no gradient to test. Every
#' derivative-free method reports a non-negative scalar of its own that goes
#' to zero as it converges, and this rule tests that. It is the default for
#' [nelder_mead()], [compass()], [sa()] and [bundle()], and a
#' gradient-based method is free to use it too.
#'
#' @param tol Numeric tolerance, a single positive number. Defaults to
#'   `1e-8`. What the tolerance is *on* depends on the method; see below.
#'
#' @details
#' A gradient-based method detects its arrival through the vanishing of
#' \eqn{\nabla f}.
#' None of the derivative-free methods can use that test, and for the
#' non-smooth problems they exist to solve it would not be the right test even
#' if they could: at the minimum of \eqn{\lvert x \rvert} any evaluated
#' subgradient is \eqn{\pm 1}, so [crit_grad()] would never fire at
#' the solution itself.
#'
#' Each such method therefore reports a non-negative scalar of its own that goes
#' to zero as it converges, and this rule tests that. What the scalar
#' *is* differs, deliberately, because the natural measure differs:
#' \describe{
#'   \item{[nelder_mead()]}{the diameter of the simplex, so the
#'     tolerance is on the parameter scale.}
#'   \item{[compass()]}{the poll size \eqn{\Delta}. This is the
#'     rule with a theorem behind it: the limit points of a pattern search with
#'     \eqn{\Delta \to 0} are Clarke stationary.}
#'   \item{[bundle()]}{the optimality estimate
#'     \eqn{\lVert p \rVert^2 + \alpha}, which vanishes exactly when zero lies
#'     in the convex hull of the collected subgradients with no linearization
#'     error. Note that this is *not* the predicted decrease, which
#'     carries a factor of the trust parameter and can therefore be driven to
#'     zero by that parameter shrinking rather than by the point becoming
#'     stationary.}
#' }
#' The measure appears in the trace as the `stationarity` column, so a run
#' can be read afterwards without knowing which method produced it.
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
#' # What the number means differs by method: a simplex diameter here, a poll
#' # size there, and the trace records whichever it is.
#' q <- function(p) sum((p - c(1, 2))^2)
#' names(minimize(nelder_mead(keep_trace = TRUE), q, c(0, 0))@trace)
#' names(minimize(compass(keep_trace = TRUE), q, c(0, 0))@trace)
#'
#' @seealso [nelder_mead()], [compass()], [sa()] and [bundle()] for the four
#'   measures, [crit_grad()] for the gradient-based rule.
#' @export
crit_stationary <- function(tol = 1e-8) {
  check_tol(tol)
  CritStationary(label = paste0("stationarity < ", format(tol)), tol = tol)
}


# --- run the budget ---------------------------------------------------------

#' @title S7 Class for the Empty Criterion
#'
#' @description
#' The rule [crit_never()] builds. Its `crit_met()` returns `FALSE` at every
#' state without reading anything, so a run carrying it ends only when a
#' budget runs out and reports `converged = FALSE`. It declares nothing
#' through [crit_needs()], so no optimizer refuses it.
#'
#' @details
#' It carries no `tol`, having nothing to compare, and is the only criterion
#' class with no property beyond `label`. [adam()] is the one shipped method
#' that defaults to it.
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

#' @title Never Stop Early
#'
#' @description
#' Builds the rule that never fires, so a run ends only when it exhausts a
#' budget and reports `converged = FALSE`. This is [adam()]'s default and is
#' the honest arrangement for a stochastic method, where every quantity a
#' convergence rule could read is an estimate.
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
#' A run that ends this way reports `converged = FALSE`, which is the
#' truth: the budget ran out, and nothing checked whether the answer was any
#' good. It is the same discipline everywhere else in the package — convergence
#' is what a rule confirmed, never what the run merely stopped doing.
#'
#' @return An S7 object of class [CritNever], inheriting from [criterion()].
#'
#' @examples
#' crit_never()
#'
#' # It fires at no state at all, however good.
#' crit_met(crit_never(), list(f_new = 0, f_old = 0, x_new = 1, x_old = 1,
#'                             gradient = c(0, 0)))
#'
#' # A run carrying it ends on its budget and says so.
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
#' The class [crit_any()] and [crit_all()] both build, and the two methods it
#' implements. Its `crit_met()` evaluates every rule it holds and reduces
#' with `any()` or `all()` according to `how`; its `crit_needs()` is the
#' union of what they need, so a combination containing a gradient rule is
#' refused by a derivative-free method exactly as the bare rule would be.
#'
#' @details
#' A `CritCombine` is itself a [criterion()], so combinations nest and the
#' label nests with them: `crit_any(crit_all(a, b), c)` reads
#' `a and b or c`. Every rule is evaluated at every call, `any()` and `all()`
#' taking the whole vector rather than short-circuiting, which costs nothing
#' worth counting against an objective evaluation.
#'
#' @param criteria A list of [criterion()] objects.
#' @param how Either `"any"` or `"all"`.
#'
#' @return An S7 object of class `CritCombine` inheriting from [criterion()],
#'   carrying `label`, `criteria` and `how`.
#'
#' @examples
#' # The needs are the union, so this is refused by a simplex method.
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
#' [crit_any()] and [crit_all()] reject the same nonsense in the same words:
#' an empty call and an argument that is not a criterion.
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
#' Combines criteria disjunctively. This is the usual arrangement: a run should
#' end as soon as any reasonable rule is satisfied.
#'
#' @details
#' **The default rule of the gradient methods.** [gd()],
#' [cg()], [bb()], [newton()],
#' [bfgs()] and [lbfgs()] default to
#' `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`: the point is
#' stationary, or the objective has stopped moving, or the parameters have.
#' Since a disjunction can only get weaker as terms are added, a run that ends
#' under this rule would have ended under a gradient rule alone at best later
#' and never earlier.
#'
#' What that buys and what it costs was measured over the package's own
#' [test_problems()], six methods on eight problems. Against the
#' gradient rule alone it converges on 44 of the 48 runs rather than 41, and
#' costs 19370 objective evaluations rather than 22299. The three it gains are
#' `cg` and `bb` on the non-smooth `abs_sum` and `gd` on
#' Beale, and none is lost.
#'
#' The cost is that a run stops sooner, so the point it reports is further from
#' the solution. Measured, 10 of the 48 end at a gradient more than a hundred
#' times larger, and the worst are runs that were reaching absurd precision
#' anyway: `bfgs` on Rosenbrock ends at 8.1e-06 rather than 4.4e-10, with
#' the objective 3.2e-13 above its minimum rather than 1.2e-21. On one run the
#' difference is real rather than cosmetic -- `cg` on `abs_sum`,
#' where the objective ends 4.5e-02 above the minimum rather than 1.9e-03, and
#' the flag reads `TRUE` where it used to read `FALSE`. That is a
#' smooth method on a non-smooth problem, where the objective stalls far from
#' the solution and a rule that reads a stall cannot tell the two apart. A
#' caller who needs stationarity asks for it: `criterion = crit_grad()`.
#'
#' [crit_rel_obj()] was in that default until 0.6.0 and is not any
#' more, because it never fired: measured over the same 48 runs, the rule with
#' it and the rule without it agree on every count, every evaluation and every
#' reported point. It remains available and useful where an objective's scale
#' is not known in advance.
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
#' Combines criteria conjunctively, for a run that should not stop until
#' several independent things agree. A conjunction can only get **stronger**
#' as terms are added, so a run ending under one would have ended under any
#' of its terms alone at the same time or earlier.
#'
#' @details
#' This is the rarer of the two combinators and is worth reaching for when
#' one rule alone is known to fire early: `crit_all(crit_grad(),
#' crit_abs_par())` asks for a stationary point at which the iterate has also
#' settled, which a run circling a flat optimum will not satisfy.
#'
#' Every rule it holds is evaluated at every iteration, so a conjunction
#' containing a rule an optimizer cannot evaluate is refused exactly as the
#' bare rule would be.
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
#' Shows a criterion in one line, as `<criterion>` and its label. The label
#' is what the result reports in `criterion_met`, so printing a rule shows
#' exactly the string a converged run will carry.
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
#' # And it is the string the result reports.
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
#' constructor rejects the same nonsense in the same words. Zero is refused
#' as well as a negative: a rule with a tolerance of zero can never fire.
#'
#' @details
#' It is also called from outside this file, by [bundle()] for its proximity
#' weights and by [bb()] for its step bounds, so the message names `tol`
#' rather than the caller's own argument. A reader who wrote
#' `bb(alpha0 = 0)` sees a message about `tol`.
#'
#' @param tol The value supplied.
#'
#' @return Invisibly `TRUE`. Raises an error naming `tol` otherwise.
#'
#' @keywords internal
check_tol <- function(tol) {
  if (length(tol) != 1L || !is.numeric(tol) || is.na(tol) || tol <= 0) {
    stop("'tol' must be a single positive number.", call. = FALSE)
  }
  invisible(TRUE)
}
