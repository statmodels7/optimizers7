#' @include criterion.R
NULL

#' @title S7 Class for Line Searches
#'
#' @description
#' This class is the abstract parent of the three line searches. A line
#' search determines how far to travel along a direction that the method has
#' already computed. It is an object, so the descent methods ([newton()], [bfgs()], [lbfgs()],
#' [cg()], [gd()], [bb()] and [chebyshev()]) share one implementation, and a
#' caller can replace it without changing the method.
#'
#' @section Notation:
#' \eqn{x} is the current point, \eqn{g = \nabla f(x)} the gradient there,
#' \eqn{d} a direction satisfying \eqn{g^\top d < 0}, and \eqn{\alpha > 0}
#' the step length that the search returns. The vector \eqn{\alpha d} is the
#' step taken; \eqn{s} elsewhere in this package is the secant vector
#' \eqn{x^{+} - x} of a quasi-Newton update, a different quantity.
#'
#' @details
#' Every subclass guarantees *sufficient decrease*. For [armijo()] and
#' [wolfe()] this is Armijo's condition
#'
#' \deqn{f(x + \alpha d) \le f(x) + c_1 \alpha\, g^\top d,
#'   \qquad 0 < c_1 < 1,}
#'
#' which requires at least the fraction \eqn{c_1} of the decrease that the
#' linear model predicts, and so rules out steps that reduce the objective by
#' an amount that vanishes faster than the step itself. [wolfe()] adds the
#' strong curvature condition
#'
#' \deqn{\lvert \nabla f(x + \alpha d)^\top d \rvert \le
#'       c_2 \lvert g^\top d \rvert, \qquad c_1 < c_2 < 1,}
#'
#' which excludes steps too short to have changed the directional derivative,
#' and which gives [bfgs()] a secant pair with usable curvature.
#' [nonmonotone()] keeps the form of Armijo's condition and replaces
#' \eqn{f(x)} by the largest of the last few values.
#'
#' The class is abstract; use [armijo()], [wolfe()] or [nonmonotone()]. A
#' fourth search would need a branch in the compiled loop, so, unlike the
#' criteria, the set of line searches cannot be extended from outside the
#' package.
#'
#' @param label A short character label, shown when the optimizer that
#'   contains the search is printed.
#'
#' @return An S7 object of class `line_search`, carrying `label`. The class
#'   is abstract, so every value is an object of one of its three subclasses.
#'
#' @examples
#' # The class is abstract: use one of the constructors.
#' try(line_search(label = "mine"))
#'
#' armijo()
#' wolfe()
#'
#' nonmonotone()
#'
#' # A method uses the search that it is given. On a smooth problem the two
#' # cost about the same; they differ on a method that needs the curvature
#' # condition and on a non-smooth problem.
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#' rbind(armijo = unlist(minimize(bfgs(line_search = armijo()), f,
#'                                c(-1.2, 1), gr = g)@counts),
#'       wolfe  = unlist(minimize(bfgs(line_search = wolfe()), f,
#'                                c(-1.2, 1), gr = g)@counts))
#'
#' # On the non-smooth problem of the battery, BFGS converges with Armijo
#' # and stops without converging with Wolfe.
#' k <- test_problems("abs_sum")[[1]]
#' rbind(armijo = c(value = minimize(bfgs(line_search = armijo()), k$fn,
#'                                   k$par, gr = k$gr)@value),
#'       wolfe  = c(value = minimize(bfgs(line_search = wolfe()), k$fn,
#'                                   k$par, gr = k$gr)@value))
#'
#' @seealso [armijo()], [wolfe()], [nonmonotone()]
#' @export
line_search <- S7::new_class(
  "line_search",
  properties = list(label = S7::class_character),
  abstract = TRUE
)


#' @title S7 Class for Armijo Backtracking
#'
#' @description
#' This class is a line search that starts from the full step and multiplies
#' it by `shrink` until the objective decreases by enough. [armijo()] builds
#' it. The search evaluates the objective at trial points and never the
#' gradient, which makes it the cheaper search per iteration.
#'
#' @details
#' Beyond the `label` of every line search, an `ArmijoSearch` holds `c1`,
#' `shrink`, `max_step` and `resolution`. It has no `c2`, since it has no
#' curvature condition; [line_search_spec()] fills that field with a value
#' that the compiled code ignores, so that the three searches are passed to
#' the compiled code as lists of the same shape.
#'
#' @param c1 The sufficient-decrease constant.
#' @param shrink The factor by which the step is multiplied on each
#'   backtrack.
#' @param max_step The most trial points allowed, including the full step.
#'   It is a count, not a length.
#' @param resolution The smallest difference in the objective that the
#'   objective can resolve, a number or a function of no arguments.
#'
#' @return An S7 object of class `ArmijoSearch` inheriting from
#'   [line_search()], with the four properties above beside `label`.
#'
#' @seealso [armijo()] for the constructor, [WolfeSearch] and
#'   [NonmonotoneSearch] for the other two.
#' @name ArmijoSearch-class
#' @aliases ArmijoSearch
#' @keywords internal
ArmijoSearch <- S7::new_class("ArmijoSearch", parent = line_search,
  properties = list(c1 = S7::class_numeric, shrink = S7::class_numeric,
                    max_step = S7::class_numeric,
                    resolution = S7::class_any))


#' @title Backtracking Line Search with the Armijo Condition
#'
#' @description
#' Multiplies the step by `shrink` until the objective decreases by enough:
#' \deqn{f(x + \alpha d) \leq f(x) + c_1 \alpha\, g^\top d .}
#' It evaluates the objective at each trial point and never the gradient,
#' which makes it the cheaper line search, suited to a method that keeps no
#' curvature approximation.
#'
#' @param c1 Sufficient-decrease constant, strictly inside \eqn{(0, 1)}.
#'   Defaults to `1e-4`, the conventional value: it requires a decrease, but
#'   only a small fraction of the decrease that the linear model predicts.
#' @param shrink Factor applied on each backtrack, strictly inside
#'   \eqn{(0, 1)}. Defaults to `0.5`.
#' @param max_step Maximum number of trial points, including the full step,
#'   before the search returns without a step; a positive whole number.
#'   Defaults to 30. The name suggests a length, but the quantity is a count.
#' @param resolution The smallest difference in the objective that the
#'   objective can resolve, in the units of the objective, or a function of
#'   no arguments that returns it when it changes during the run. Defaults to
#'   `0`, which disables the test. See the section below.
#'
#' @section What the objective can resolve:
#' An objective computed by a procedure, and not by a formula, can return
#' slightly different values for the same argument: a fit warm-started from
#' the point where the last evaluation ended, a quadrature whose panels move,
#' a simulation. Differences below that spread carry no information, and a
#' search that has to verify a smaller decrease backtracks until it has used
#' `max_step` trial points.
#'
#' The test is made once, before the first trial point, on the decrease
#' that the linear model of the method predicts over the full step,
#' \eqn{\alpha_0 \lvert g^\top d\rvert}. When that decrease is below
#' `resolution` the search returns at once, and the run reports that the
#' point is optimal to the accuracy of the objective. That is a weaker
#' statement than the one made when a stopping rule is met, and the result
#' reports it with a different label.
#'
#' The test is not applied inside the backtracking loop. There the two
#' situations cannot be distinguished, since \eqn{x + \alpha d \to x} as the
#' step shrinks and the objective stops resolving the change whether the
#' point is optimal or the direction is wrong. At the full step they differ:
#' a wrong direction predicts a large decrease, does not obtain it, and the
#' search reports a failure. Inside the loop the test would report an
#' incorrect gradient at a point far from stationarity as convergence.
#'
#' The quantity is the predicted decrease and not the Armijo threshold
#' \eqn{c_1 \alpha_0 \lvert g^\top d\rvert}, which at the default
#' \eqn{c_1 = 10^{-4}} is four orders of magnitude smaller and would fire
#' while the method can still make progress.
#'
#' ## A resolution that changes during the run
#'
#' When the objective becomes more precise as the run proceeds (for example
#' a fit warm-started from the previous evaluation), the resolution measured
#' at the start is the worst value of the whole run. A function of no
#' arguments passed instead of a number is called at every iteration, once
#' per invocation of the search and not once per trial, so it costs one call
#' per iteration. Its value is the resolution for the step about to be
#' taken; a value that is not finite and positive is read as `0`, which
#' disables the test.
#'
#' @details
#' The term \eqn{c_1 \alpha\, g^\top d} is required. A test that requires
#' only \eqn{f_{new} \le f} accepts steps that leave the objective
#' unchanged: on \eqn{f(x) = \lVert x \rVert^2}, whose Hessian is \eqn{2I},
#' a unit step along \eqn{-g} reflects the iterate through the minimum and
#' leaves the objective *exactly* unchanged. The weaker test accepts that
#' step, the iterate oscillates indefinitely, and a stopping rule that reads
#' the objective sees no change and reports convergence at a point that is
#' not a minimum.
#'
#' A quasi-Newton method also needs the curvature condition, which
#' [wolfe()] enforces.
#'
#' @return An S7 object of class [ArmijoSearch], inheriting from
#'   [line_search()], to be passed as an optimizer's `line_search`.
#'
#' @examples
#' armijo()
#' armijo(shrink = 0.2, max_step = 10)
#'
#' minimize(gd(line_search = armijo(shrink = 0.2)),
#'          function(p) sum((p - 1:2)^2), c(0, 0))@par
#'
#' # Barzilai-Borwein relies on steps that go uphill. An Armijo condition
#' # rejects every such step, and the run costs more than twice as many
#' # evaluations as under the default nonmonotone search.
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#' c(armijo      = minimize(bb(line_search = armijo()), f, c(-1.2, 1),
#'                          gr = g)@counts[["f"]],
#'   nonmonotone = minimize(bb(), f, c(-1.2, 1), gr = g)@counts[["f"]])
#'
#' # Constants outside their intervals are rejected, and the error names
#' # the argument.
#' try(armijo(c1 = 1))
#' try(armijo(shrink = 1.5))
#' try(armijo(max_step = 2.5))
#'
#' @seealso [wolfe()] for the curvature condition, [nonmonotone()] for the
#'   relaxed reference value, [line_search()] for the parent class.
#' @references
#' Armijo, L. (1966). Minimization of functions having Lipschitz
#' continuous first partial derivatives. *Pacific Journal of
#' Mathematics* **16**, 1--3.
#'
#' Zoutendijk, G. (1970). Nonlinear programming, computational methods.
#' In J. Abadie (ed.), *Integer and Nonlinear Programming*,
#' 37--86. North-Holland, Amsterdam.
#'
#' @export
armijo <- function(c1 = 1e-4, shrink = 0.5, max_step = 30, resolution = 0) {
  check_unit(c1, "c1")
  check_unit(shrink, "shrink")
  check_count(max_step, "max_step")
  check_resolution(resolution)
  ArmijoSearch(label = paste0("Armijo backtracking (c1 = ", format(c1), ")"),
               c1 = c1, shrink = shrink, max_step = max_step,
               resolution = resolution)
}


#' Check a Line Search's Resolution
#'
#' @description
#' Checks that `resolution` is a single non-negative finite number, where
#' zero disables the test, or a function of no arguments that returns one.
#'
#' @details
#' The function form is for an objective whose resolution changes during the
#' run. It is called once per invocation of the search and not once per
#' trial, so it costs one call per iteration. A caller whose objective
#' becomes more precise during the run can then return the resolution at the
#' current point instead of the value measured at the start.
#'
#' @param x What the constructor was given.
#'
#' @return `NULL`, invisibly; called for the error.
#'
#' @seealso [armijo()]
#'
#' @keywords internal
check_resolution <- function(x) {
  if (is.function(x)) {
    if (length(formals(x))) {
      stop(paste0("a 'resolution' given as a function must take no arguments;",
                  " it is called\n  for the resolution in force at the",
                  " iteration about to be taken."), call. = FALSE)
    }
    return(invisible(NULL))
  }
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || x < 0 ||
      !is.finite(x)) {
    stop(paste0("'resolution' must be a single non-negative finite number,",
                " 0 to disable the test,\n  or a function of no arguments",
                " returning one."), call. = FALSE)
  }
  invisible(NULL)
}


#' @title S7 Class for the Strong Wolfe Line Search
#'
#' @description
#' This class is a line search that brackets a step satisfying both Wolfe
#' conditions and then bisects inside the bracket. [wolfe()] builds it. The
#' search evaluates the gradient at trial points as well as the objective,
#' which is the cost of the curvature condition that a quasi-Newton method
#' needs.
#'
#' @details
#' Beyond the `label` of every line search, a `WolfeSearch` holds `c1`,
#' `c2`, `max_step` and `resolution`. It has no `shrink`, since it takes no
#' backtracking steps of a fixed ratio; [line_search_spec()] fills that field
#' with a value that the compiled code ignores.
#'
#' @param c1 The sufficient-decrease constant.
#' @param c2 The curvature constant.
#' @param max_step The most trial steps allowed in each phase. It is a
#'   count.
#' @param resolution The smallest difference in the objective that the
#'   objective can resolve, a number or a function of no arguments.
#'
#' @return An S7 object of class `WolfeSearch` inheriting from
#'   [line_search()], with the four properties above beside `label`.
#'
#' @seealso [wolfe()] for the constructor, [ArmijoSearch] and
#'   [NonmonotoneSearch] for the other two.
#' @name WolfeSearch-class
#' @aliases WolfeSearch
#' @keywords internal
WolfeSearch <- S7::new_class("WolfeSearch", parent = line_search,
  properties = list(c1 = S7::class_numeric, c2 = S7::class_numeric,
                    max_step = S7::class_numeric,
                    resolution = S7::class_any))


#' @title Line Search Satisfying the Strong Wolfe Conditions
#'
#' @description
#' Searches for a step that satisfies both Wolfe conditions, sufficient
#' decrease and curvature:
#' \deqn{f(x + \alpha d) \leq f(x) + c_1 \alpha\, g^\top d, \qquad
#'       \lvert g(x + \alpha d)^\top d \rvert \leq c_2 \lvert g^\top d \rvert .}
#' It evaluates the gradient at trial points as well as the objective, which
#' is the cost of the second condition, the one that [bfgs()] and [lbfgs()]
#' need.
#'
#' @param c1 Sufficient-decrease constant, strictly inside \eqn{(0, 1)}.
#'   Defaults to `1e-4`.
#' @param c2 Curvature constant, with \eqn{c_1 < c_2 < 1}. Defaults to `0.9`,
#'   the usual choice for a quasi-Newton method; [cg()] uses `0.1`, a more
#'   exact line search. A `c2` at or below `c1` is rejected, because a step
#'   that satisfies both conditions is guaranteed to exist only when
#'   \eqn{c_1 < c_2 < 1}.
#' @param max_step Maximum trial steps in each of the bracketing and zoom
#'   phases, a positive whole number. Defaults to 30.
#' @param resolution The smallest difference in the objective that the
#'   objective can resolve, in the units of the objective, or a function of
#'   no arguments that returns it when it changes during the run. Defaults to
#'   `0`, which disables the test; see [armijo()] for the purpose of the test
#'   and for the reason why it is applied at the full step, before the search
#'   begins.
#'
#' @details
#' # The curvature condition
#'
#' [armijo()] does not enforce it, and a quasi-Newton method relies on it.
#' BFGS builds its approximation from the secant pair \eqn{(s, y)} with
#' \eqn{s = \alpha d} and \eqn{y = g_{new} - g_{old}}, and a step so short
#' that the gradient has barely changed gives a pair with no usable
#' curvature: the update is then skipped, or it degrades the approximation.
#' Requiring the directional derivative to have shrunk by the factor
#' \eqn{c_2} rules this out. On the problems of [test_problems()], BFGS
#' under [armijo()] skips an update on a few problems only; see [bfgs()].
#'
#' # The implementation
#'
#' The search brackets a step, doubling the trial step until the bracket is
#' found, and then bisects inside the bracket instead of interpolating a
#' polynomial. Bisection costs a few more evaluations and does not depend on
#' the shape of the function inside the bracket. Each phase is limited to
#' `max_step` trials, and a search that exhausts them returns without a
#' step. When the bracket collapses before the curvature condition holds,
#' the end of the bracket that satisfies sufficient decrease is returned and
#' the trace records `curvature condition not met`.
#'
#' The gradient evaluations at the trial points make this search more costly
#' per iteration than [armijo()]; the total cost of a run depends on the
#' problem.
#'
#' @return An S7 object of class [WolfeSearch], inheriting from
#'   [line_search()], to be passed as an optimizer's `line_search`.
#'
#' @examples
#' wolfe()
#' wolfe(c2 = 0.1)
#'
#' # Gradient descent with this search comes close to the minimum (1, 1) in
#' # 2000 iterations.
#' minimize(gd(line_search = wolfe(), maxit = 2000),
#'          function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2,
#'          c(-1.2, 1))@par
#'
#' # The smaller constant matters for conjugate gradients, which has no
#' # curvature approximation to compensate for an inexact step.
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#' c(tight = minimize(cg(line_search = wolfe(c2 = 0.1)), f, c(-1.2, 1),
#'                    gr = g)@iterations,
#'   loose = minimize(cg(line_search = wolfe(c2 = 0.9)), f, c(-1.2, 1),
#'                    gr = g)@iterations)
#'
#' # c2 must exceed c1.
#' try(wolfe(c1 = 0.5, c2 = 0.4))
#'
#' @seealso [armijo()] for the cheap search, [nonmonotone()] for the relaxed
#'   reference value, [bfgs()] for the method that needs this one.
#' @references
#' Wolfe, P. (1969). Convergence conditions for ascent methods.
#' *SIAM Review* **11**, 226--235.
#'
#' Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*,
#' 2nd edition. Springer, New York.
#'
#' @export
wolfe <- function(c1 = 1e-4, c2 = 0.9, max_step = 30, resolution = 0) {
  check_unit(c1, "c1")
  check_unit(c2, "c2")
  if (c2 <= c1) {
    stop("'c2' must be greater than 'c1'; only then is a step that ",
         "satisfies both conditions guaranteed to exist.", call. = FALSE)
  }
  check_count(max_step, "max_step")
  check_resolution(resolution)
  WolfeSearch(label = paste0("strong Wolfe (c1 = ", format(c1),
                             ", c2 = ", format(c2), ")"),
              c1 = c1, c2 = c2, max_step = max_step,
              resolution = resolution)
}


#' Describe a Line Search to the C++ Side
#'
#' @description
#' Converts a [line_search()] object into the plain list that the compiled
#' loop reads, so the C++ code does not have to handle S7 objects. Every
#' search produces the same seven fields, whatever its subclass. The compiled
#' search branches on `type`, and the descent loop reads `memory` to form the
#' reference value.
#'
#' @details
#' A field that the search does not use is filled with a value that the
#' compiled code ignores: `c2 = 0.9` for the two backtracking searches,
#' `shrink = 0.5` for Wolfe, `memory = 0` for the two monotone ones. The
#' shape is fixed, so one C++ structure reads all three.
#'
#' There are only two types. [nonmonotone()] is passed as
#' `type = "armijo"` with its `memory`, since it differs from [armijo()] only
#' in the reference value, so the compiled search does not need a third
#' branch.
#' With `memory = 0` the two lists are identical.
#'
#' @param x A [line_search()] object.
#'
#' @return A list of seven: `type` (`"armijo"` or `"wolfe"`), `c1`, `c2`,
#'   `shrink`, `max_step` (integer), `memory` (integer) and `resolution`.
#'
#' @aliases line_search_spec.ArmijoSearch line_search_spec.WolfeSearch
#' @keywords internal
line_search_spec <- S7::new_generic("line_search_spec", "x",
                                    function(x) S7::S7_dispatch())

S7::method(line_search_spec, ArmijoSearch) <- function(x) {
  list(type = "armijo", c1 = x@c1, c2 = 0.9, shrink = x@shrink,
       max_step = as.integer(x@max_step), memory = 0L,
       resolution = x@resolution)
}

S7::method(line_search_spec, WolfeSearch) <- function(x) {
  list(type = "wolfe", c1 = x@c1, c2 = x@c2, shrink = 0.5,
       max_step = as.integer(x@max_step), memory = 0L,
       resolution = x@resolution)
}


#' @title S7 Class for the Nonmonotone Line Search
#'
#' @description
#' This class is Armijo backtracking whose reference value is the largest of
#' the last `memory + 1` objective values, including the current one, so a
#' step may increase the objective. [nonmonotone()] builds it.
#'
#' @details
#' Beyond the `label` of every line search, a `NonmonotoneSearch` holds
#' `c1`, `shrink`, `memory`, `max_step` and `resolution`. Only this class has
#' the property `memory`, and at `memory = 0` the object behaves exactly as
#' an [ArmijoSearch] with the same constants.
#'
#' @param c1 The sufficient-decrease constant.
#' @param shrink The factor by which the step is multiplied on each
#'   backtrack.
#' @param memory The number of earlier values, besides the current one, over
#'   which the reference is taken.
#' @param max_step The most trial points allowed, including the full step.
#'   It is a count.
#' @param resolution The smallest difference in the objective that the
#'   objective can resolve, a number or a function of no arguments.
#'
#' @return An S7 object of class `NonmonotoneSearch` inheriting from
#'   [line_search()], with the five properties above beside `label`.
#'
#' @seealso [nonmonotone()] for the constructor, [bb()] for the method that
#'   uses it by default.
#' @name NonmonotoneSearch-class
#' @aliases NonmonotoneSearch
#' @keywords internal
NonmonotoneSearch <- S7::new_class("NonmonotoneSearch", parent = line_search,
  properties = list(
    c1         = S7::class_numeric,
    shrink     = S7::class_numeric,
    memory     = S7::class_numeric,
    max_step   = S7::class_numeric,
    resolution = S7::class_any
  ))


#' @title Nonmonotone Backtracking
#'
#' @description
#' Performs Armijo backtracking that compares the trial value with the
#' largest of the last few objective values and not with the current one, so
#' a step may increase the objective.
#'
#' @param c1 Sufficient-decrease constant. Defaults to `1e-4`.
#' @param shrink Factor applied to the step on each backtrack. Defaults to
#'   `0.5`.
#' @param memory The number of earlier values, besides the current one, over
#'   which the reference is taken. Defaults to `10`; `0` makes the search
#'   identical to [armijo()].
#' @param max_step Maximum number of trial points, including the full step,
#'   before the search returns without a step. Defaults to `30`.
#' @param resolution The smallest difference in the objective that the
#'   objective can resolve, in the units of the objective, or a function of
#'   no arguments that returns it when it changes during the run. Defaults to
#'   `0`, which disables the test; see [armijo()] for the purpose of the test
#'   and for the reason why it is applied at the full step, before the search
#'   begins.
#'
#' @details
#' The condition is that of Grippo, Lampariello and Lucidi:
#' \deqn{f(x_k + \alpha d_k) \le \max_{0 \le j \le m} f(x_{k-j}) +
#'   c_1 \alpha\, g_k^\top d_k,}
#' which is Armijo's condition with the reference replaced by the largest of
#' the last \eqn{m+1} values. An accepted step improves on the largest of the
#' recent values, but it need not improve on the current one.
#'
#' # Purpose
#'
#' Some methods are efficient *because* of steps that increase the
#' objective. [bb()] is the main case: its step length is a curvature
#' estimate taken from the last secant pair, and following that estimate
#' sometimes leads to a higher value. An Armijo condition rejects every step
#' that increases the objective and backtracks until it finds a shorter one,
#' which removes the steps on which the method relies. On the Rosenbrock
#' function [bb()] needs fewer than half as many evaluations of the objective
#' under `nonmonotone()` as under [armijo()] (see the examples).
#'
#' The cost is a weaker guarantee. A monotone method cannot cycle, since the
#' objective is a decreasing sequence bounded below; a nonmonotone method
#' relies on the finite memory for that, and its convergence result is more
#' delicate. The search is intended for methods such as [bb()], whose steps
#' are not required to decrease the objective.
#'
#' # memory = 0
#'
#' At `memory = 0` the reference is the current value and the condition is
#' exactly Armijo's, so `nonmonotone(memory = 0)` and [armijo()] produce the
#' same run, and a comparison between the two settings isolates the effect
#' of the memory.
#'
#' # No nonmonotone Wolfe search
#'
#' The curvature condition concerns the gradient at the trial point and does
#' not involve the value against which the decrease is measured. A
#' nonmonotone Wolfe search would therefore be a fourth object and not an
#' option of this one, and the package does not provide one.
#'
#' @return A [line_search()] object.
#'
#' @references
#' Grippo, L., Lampariello, F. and Lucidi, S. (1986). A nonmonotone line search
#' technique for Newton's method. *SIAM Journal on Numerical Analysis*
#' **23**, 707--716.
#'
#' Raydan, M. (1997). The Barzilai and Borwein gradient method for the large
#' scale unconstrained minimization problem. *SIAM Journal on
#' Optimization* **7**, 26--33.
#'
#' @examples
#' nonmonotone()
#' nonmonotone(memory = 5)
#'
#' # The example shows the effect on bb(), the method that uses this search
#' # by default.
#' f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
#' gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
#'                     200 * (p[2] - p[1]^2))
#' evals <- function(ls) minimize(bb(line_search = ls), f, c(-1.2, 1),
#'                                gr = gr)@counts[["f"]]
#' c(nonmonotone = evals(nonmonotone()), armijo = evals(armijo()))
#'
#' # With memory = 0 the search is armijo(), so the comparison above
#' # isolates the memory.
#' identical(evals(nonmonotone(memory = 0)), evals(armijo()))
#'
#' @seealso [armijo()] for the monotone version, [bb()] for the method that
#'   uses this one by default, [line_search()] for the parent class.
#' @export
nonmonotone <- function(c1 = 1e-4, shrink = 0.5, memory = 10, max_step = 30,
                        resolution = 0) {
  check_unit(c1, "c1")
  check_unit(shrink, "shrink")
  check_count(max_step, "max_step")
  check_resolution(resolution)
  check_whole(memory, "memory")
  NonmonotoneSearch(
    label = paste0("nonmonotone backtracking (memory = ", format(memory), ")"),
    c1 = c1, shrink = shrink, memory = memory, max_step = max_step,
    resolution = resolution)
}


#' @rdname line_search_spec
#' @name line_search_spec.NonmonotoneSearch
#' @keywords internal
S7::method(line_search_spec, NonmonotoneSearch) <- function(x) {
  list(type = "armijo", c1 = x@c1, c2 = 0.9, shrink = x@shrink,
       max_step = as.integer(x@max_step), memory = as.integer(x@memory),
       resolution = x@resolution)
}

#' @title Print Method for Line Searches
#' @name print.line_search
#' @param x A [line_search()] object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' print(wolfe())
#' @keywords internal
S7::method(print, line_search) <- function(x, ...) {
  cat("<line_search> ", x@label, "\n", sep = "")
  invisible(x)
}


#' Validate a Constant in the Unit Interval
#'
#' @description
#' Checks that the value is a single number strictly inside
#' \eqn{(0, 1)}. Both endpoints are rejected: `c1 = 0` requires no decrease
#' at all, and `c1 = 1` requires the whole decrease that the linear model
#' predicts, which a strictly convex objective cannot supply.
#'
#' @param v The value.
#' @param nm Its name, for the message.
#'
#' @return Invisibly `TRUE`. Raises an error naming `nm` otherwise.
#'
#' @keywords internal
check_unit <- function(v, nm) {
  if (length(v) != 1L || !is.numeric(v) || is.na(v) || v <= 0 || v >= 1) {
    stop("'", nm, "' must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a Positive Whole Number
#'
#' @description
#' Checks that the value is a single positive whole number no larger than
#' `.Machine$integer.max`. It is used for `max_step`, which counts trial
#' points, and for the `memory` of [lbfgs()], which counts secant pairs. The
#' count is passed to the compiled code as an integer, so a fractional value is
#' rejected (`memory = 2.5` would otherwise run at a memory of 2), and so is
#' an infinite one, which `as.integer()` would turn into `NA`.
#'
#' @param v The value.
#' @param nm Its name, for the message.
#'
#' @return Invisibly `TRUE`. Raises an error naming `nm` otherwise.
#'
#' @seealso [check_whole()] for a count that may be zero.
#' @keywords internal
check_count <- function(v, nm) {
  if (length(v) != 1L || !is.numeric(v) || is.na(v) || v < 1 ||
      v > .Machine$integer.max || v != round(v)) {
    stop("'", nm, "' must be a single positive whole number.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a Non-Negative Whole Number
#'
#' @description
#' Checks that the value is a single whole number, including zero, no larger
#' than `.Machine$integer.max`. It is used where a count gives how many of
#' something are tolerated, so that zero means that every one is rejected:
#' the `memory` of [nonmonotone()], where zero makes the reference the
#' current value and the search identical to [armijo()], and the `max_skip`
#' of [bfgs()], where zero resets the approximation on the first skipped
#' update.
#'
#' A negative value is rejected because the compiled code treats it as zero:
#' it compares `max_skip` with `>=` and tests `memory > 0`. An infinite value
#' is rejected because `as.integer()` would turn it into `NA`.
#'
#' @param v The value.
#' @param nm Its name, for the message.
#'
#' @return Invisibly `TRUE`. Raises an error naming `nm` otherwise.
#'
#' @seealso [check_count()] for a count that has to be at least one.
#' @keywords internal
check_whole <- function(v, nm) {
  if (length(v) != 1L || !is.numeric(v) || is.na(v) || v < 0 ||
      v > .Machine$integer.max || v != round(v)) {
    stop("'", nm, "' must be a single non-negative whole number.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Validate a Non-Negative Threshold
#'
#' @description
#' Checks that the value is a single finite number that is not negative.
#' It is used for the curvature threshold `curvature_tol` of [bfgs()] and
#' [lbfgs()], which is compared as \eqn{s'y > \mathrm{curv\_tol} \, \lVert s \rVert
#' \lVert y \rVert}. Zero is admitted and is the textbook curvature
#' condition \eqn{s'y > 0}. A negative value is rejected because it admits
#' pairs with \eqn{s'y \le 0}, which weakens the protection given by the
#' skipped updates; by the Cauchy-Schwarz inequality a value below \eqn{-1}
#' admits every pair.
#'
#' @param v The value.
#' @param nm Its name, for the message.
#'
#' @return Invisibly `TRUE`. Raises an error naming `nm` otherwise.
#'
#' @keywords internal
check_nonneg <- function(v, nm) {
  if (length(v) != 1L || !is.numeric(v) || is.na(v) || !is.finite(v) ||
      v < 0) {
    stop("'", nm, "' must be a single finite number that is not negative.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Fetched rather than captured: comparing S7 classes by identity breaks when
# the class object is re-created.
#' The line_search Class Object
#'
#' @description
#' Returns the class object from the namespace at the time of the call, so
#' that a class check does not depend on whether the class object is the one
#' created when the package was loaded.
#'
#' @details
#' Comparing S7 classes is a comparison of objects, so it is `FALSE` for a
#' class rebuilt from the same definition, which happens under any loader
#' that re-evaluates the code instead of loading it, \pkg{covr} among
#' them.
#'
#' @return The [line_search()] class object.
#'
#' @keywords internal
line_search_class <- function() line_search
