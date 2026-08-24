#' @include criterion.R
NULL

#' @title S7 Class for Line Searches
#'
#' @description
#' The abstract parent of the three line searches. A line search answers one
#' question, how far to travel along a direction another piece of the method
#' has already chosen, and it answers it as an object, so Newton, BFGS,
#' L-BFGS, conjugate gradients, gradient descent and Barzilai-Borwein share
#' one carefully written answer and a caller can replace it without touching
#' the method.
#'
#' @section Notation:
#' \eqn{x} is the current point, \eqn{g = \nabla f(x)} the gradient there,
#' \eqn{d} a direction satisfying \eqn{g^\top d < 0}, and \eqn{\alpha > 0}
#' the **step length** the search returns. The vector \eqn{\alpha d} is the
#' step taken; \eqn{s} elsewhere in this package is the secant vector
#' \eqn{x^{+} - x} of a quasi-Newton update, a different quantity.
#'
#' @details
#' What every subclass guarantees is *sufficient decrease*, Armijo's
#' condition
#'
#' \deqn{f(x + \alpha d) \le f(x) + c_1 \alpha\, g^\top d,
#'   \qquad 0 < c_1 < 1,}
#'
#' which asks for a fraction \eqn{c_1} of the decrease the linear model
#' predicts, and so rules out steps that shrink the objective by an amount
#' vanishing faster than the step itself. [wolfe()] adds the strong
#' curvature condition
#'
#' \deqn{\lvert \nabla f(x + \alpha d)^\top d \rvert \le
#'       c_2 \lvert g^\top d \rvert, \qquad c_1 < c_2 < 1,}
#'
#' which excludes steps too short to have moved the directional derivative,
#' and is the guarantee [bfgs()] needs for its secant pair to carry usable
#' curvature. [nonmonotone()] keeps Armijo's condition and replaces \eqn{f(x)}
#' by the worst of the last few values.
#'
#' The class is abstract; use [armijo()], [wolfe()] or [nonmonotone()]. A
#' fourth would need a branch in the compiled loop, so the set is not
#' extensible from outside the package as the criteria are.
#'
#' @param label A short character label, shown when the optimizer carrying it
#'   is printed.
#'
#' @return An S7 object of class `line_search`, carrying `label`. The class
#'   is abstract, so every value is an object of one of its three subclasses.
#'
#' @examples
#' # Abstract: use one of the constructors.
#' try(line_search(label = "mine"))
#'
#' armijo()
#' wolfe()
#'
#' nonmonotone()
#'
#' # A method takes whichever it is given. On a problem where nothing goes
#' # wrong the two cost about the same; where they differ is on a method that
#' # needs the curvature condition, and on how a non-smooth problem is
#' # handled.
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#' rbind(armijo = unlist(minimize(bfgs(line_search = armijo()), f,
#'                                c(-1.2, 1), gr = g)@counts),
#'       wolfe  = unlist(minimize(bfgs(line_search = wolfe()), f,
#'                                c(-1.2, 1), gr = g)@counts))
#'
#' # The non-smooth problem of the battery, where Armijo is the better of the
#' # two and Wolfe stops without converging.
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
#' A line search that starts from the full step and halves it until the
#' objective falls by enough. Built by [armijo()]. It evaluates the objective
#' at trial points and never the gradient, which makes it the cheap choice
#' per iteration.
#'
#' @details
#' Beyond the `label` every line search carries, an `ArmijoSearch` holds
#' `c1`, `shrink`, `max_step` and `resolution`. It has no `c2`, having no
#' curvature condition; [line_search_spec()] fills that field with a value
#' the compiled side ignores so that all three searches describe themselves
#' in the same shape.
#'
#' @param c1 The sufficient-decrease constant.
#' @param shrink The factor the step is multiplied by on each backtrack.
#' @param max_step The most backtracks allowed. A **count**, not a length.
#' @param resolution What the objective can tell apart, a number or a
#'   function of no arguments.
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
#' Halves the step until the objective decreases by enough:
#' \deqn{f(x + \alpha d) \leq f(x) + c_1 \alpha\, g^\top d .}
#' Evaluates the objective at each trial point and never the gradient, which
#' makes it the cheap line search and the one a method with no curvature
#' approximation to protect should use.
#'
#' @param c1 Sufficient-decrease constant, strictly inside \eqn{(0, 1)}.
#'   Defaults to `1e-4`, the conventional value: it demands a decrease, but
#'   only a tiny fraction of what the linear model predicts, so it almost
#'   never rejects a sensible step.
#' @param shrink Factor applied on each backtrack, strictly inside
#'   \eqn{(0, 1)}. Defaults to `0.5`.
#' @param max_step Maximum **backtracks** before the search gives up, a
#'   positive whole number. Defaults to 30. The name says step and the
#'   quantity is a count.
#' @param resolution The smallest difference in the objective that means
#'   anything, in the objective's own units, or a function of no arguments
#'   returning it where it moves as the run goes. Defaults to `0`, which
#'   does not ask the question. See the section below.
#'
#' @section What the objective can resolve:
#' An objective computed by a procedure instead of by a formula returns
#' slightly different values for the same argument: a fit warm-started from
#' wherever the last evaluation ended, a quadrature whose panels move, a
#' simulation. Below that spread its values carry no information, and a
#' search asked to verify a smaller decrease backtracks to exhaustion.
#'
#' The test is made **once**, before any trial is paid for, on the
#' improvement the method's own linear model predicts over the full step,
#' \eqn{\alpha_0 \lvert g^\top d\rvert}. Where that is below `resolution` the
#' search returns immediately, and the run reports that the point is optimal
#' to the accuracy the objective has. That is a weaker statement than a
#' stopping rule being met, and it is reported in different words.
#'
#' **It is not asked inside the backtracking loop, and that is what makes it
#' safe.** There the two situations cannot be told apart, since
#' \eqn{x + \alpha d \to x} as the step shrinks and the objective stops
#' resolving the change whether the point is optimal or the direction is
#' wrong. Tested at the full step they separate: a bad direction predicts a
#' large improvement, does not obtain it, and is still reported as the failure
#' it is. Measured, with the test inside the loop a mis-stated gradient at a
#' point nowhere near stationary was promoted to a converged run.
#'
#' The quantity is the predicted decrease and not the Armijo demand
#' \eqn{c_1 \alpha_0 \lvert g^\top d\rvert}, which is four orders smaller and
#' would fire where the method still had real progress to make.
#'
#' # A resolution that moves
#'
#' Where the objective settles as the run goes, a fit warm-started from the
#' previous evaluation locating its own answer better each time, the
#' resolution at the start is the reading from the worst point of the whole
#' run. Passing a **function** of no arguments instead of a number has it
#' asked again at every iteration, once per invocation of the search and not
#' once per trial, so it costs one call an iteration. What the function
#' returns is the resolution in force for the step about to be taken; a value
#' that is not finite and positive is read as `0`, which asks nothing.
#'
#' @details
#' The \eqn{c_1 \alpha\, g^\top d} term is essential, and dropping it to test
#' merely \eqn{f_{new} \le f} is a real defect. On a quadratic with unit step
#' the gradient update reflects the iterate through the minimum, leaving the
#' objective *exactly* unchanged; the weak test accepts that step, the
#' iterate oscillates for ever, and a stopping rule watching the objective
#' sees no change and reports convergence at a point that is not a minimum.
#'
#' The search evaluates the objective at trial points and never the gradient.
#' That is enough for a method that has only to make progress, and not enough
#' for a quasi-Newton method, which needs [wolfe()].
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
#' # The cheapness shows where the method has to backtrack often. Under
#' # Barzilai-Borwein, which relies on steps that go uphill, an Armijo
#' # condition rejects exactly those and the run costs more than twice as
#' # much.
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#' c(armijo      = minimize(bb(line_search = armijo()), f, c(-1.2, 1),
#'                          gr = g)@counts[["f"]],
#'   nonmonotone = minimize(bb(), f, c(-1.2, 1), gr = g)@counts[["f"]])
#'
#' # Constants outside their intervals are refused by name.
#' try(armijo(c1 = 1))
#' try(armijo(shrink = 1.5))
#' try(armijo(max_step = 2.5))
#'
#' @seealso [wolfe()] for the curvature condition, [nonmonotone()] for the
#'   relaxed reference value, [line_search()] for the shared contract.
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
#' A single non-negative finite number, zero meaning the question is not asked,
#' or a function of no arguments returning one.
#'
#' @details
#' The function form is for an objective whose resolution MOVES. It is asked
#' once per invocation of the search rather than per trial, so it costs one
#' call an iteration, and it is what lets a caller whose objective settles as
#' it goes report the resolution of the current point instead of the reading
#' from the worst-located point of the run.
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
                " 0 to leave it unasked,\n  or a function of no arguments",
                " returning one."), call. = FALSE)
  }
  invisible(NULL)
}


#' @title S7 Class for the Strong Wolfe Line Search
#'
#' @description
#' A line search that brackets a step satisfying both Wolfe conditions and
#' then bisects inside the bracket. Built by [wolfe()]. It evaluates the
#' **gradient** at trial points as well as the objective, which is what the
#' curvature condition costs and what a quasi-Newton method needs.
#'
#' @details
#' Beyond the `label` every line search carries, a `WolfeSearch` holds `c1`,
#' `c2`, `max_step` and `resolution`. It has no `shrink`, taking no
#' backtracking steps of a fixed ratio; [line_search_spec()] fills that field
#' with a value the compiled side ignores.
#'
#' @param c1 The sufficient-decrease constant.
#' @param c2 The curvature constant.
#' @param max_step The most trial steps allowed in each phase. A **count**.
#' @param resolution What the objective can tell apart, a number or a
#'   function of no arguments.
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
#' Finds a step satisfying both Wolfe conditions, sufficient decrease and
#' curvature together:
#' \deqn{f(x + \alpha d) \leq f(x) + c_1 \alpha\, g^\top d, \qquad
#'       \lvert g(x + \alpha d)^\top d \rvert \leq c_2 \lvert g^\top d \rvert .}
#' Evaluates the gradient at trial points as well as the objective, which is
#' what the second condition costs and what [bfgs()] and [lbfgs()] need.
#'
#' @param c1 Sufficient-decrease constant, strictly inside \eqn{(0, 1)}.
#'   Defaults to `1e-4`.
#' @param c2 Curvature constant, with \eqn{c_1 < c_2 < 1}. Defaults to `0.9`,
#'   the usual choice for a quasi-Newton method; [cg()] defaults to `0.1`,
#'   wanting a more exact line search. A `c2` at or below `c1` is refused,
#'   the two conditions being unsatisfiable together then.
#' @param max_step Maximum trial steps in **each** of the bracketing and zoom
#'   phases, a positive whole number. Defaults to 30.
#' @param resolution The smallest difference in the objective that means
#'   anything, in the objective's own units, or a function of no arguments
#'   returning it where it moves as the run goes. Defaults to `0`, which
#'   does not ask the question; see [armijo()] for what it is for and
#'   why it is asked at the full step rather than during the search.
#'
#' @details
#' # What the curvature condition is for
#'
#' [armijo()] cannot provide it, and a quasi-Newton method needs it to work
#' at all. BFGS builds its approximation from the secant pair \eqn{(s, y)}
#' with \eqn{s = \alpha d} and \eqn{y = g_{new} - g_{old}}, and a step so
#' short that the gradient has barely moved gives a pair carrying no
#' curvature: the update is then either skipped or it corrupts the matrix.
#' Requiring the gradient along the direction to have shrunk by a factor
#' \eqn{c_2} is exactly the guarantee that this does not happen.
#'
#' In practice the guarantee is worth less than it sounds, and the honest
#' measurement belongs here: over the eight problems of [test_problems()],
#' BFGS under Armijo skips an update on three of them, once each. See
#' [bfgs()] for the numbers.
#'
#' # The implementation
#'
#' Bracketing and zoom, with **bisection** inside the zoom instead of
#' polynomial interpolation. That costs a few more evaluations and cannot be
#' defeated by an awkwardly shaped interval.
#'
#' Gradient evaluations at trial points are the price, so this is the dearer
#' choice per iteration and usually the cheaper one per problem.
#'
#' @return An S7 object of class [WolfeSearch], inheriting from
#'   [line_search()], to be passed as an optimizer's `line_search`.
#'
#' @examples
#' wolfe()
#' wolfe(c2 = 0.1)
#'
#' # Even gradient descent gets to the answer with it, given the budget.
#' minimize(gd(line_search = wolfe(), maxit = 2000),
#'          function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2,
#'          c(-1.2, 1))@par
#'
#' # The tighter constant matters to conjugate gradients, which has no
#' # curvature approximation to repair a loose step with.
#' f <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' g <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                    200 * (p[2] - p[1]^2))
#' c(tight = minimize(cg(line_search = wolfe(c2 = 0.1)), f, c(-1.2, 1),
#'                    gr = g)@iterations,
#'   loose = minimize(cg(line_search = wolfe(c2 = 0.9)), f, c(-1.2, 1),
#'                    gr = g)@iterations)
#'
#' # The two conditions have to be satisfiable together.
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
    stop("'c2' must be greater than 'c1'; the conditions are unsatisfiable ",
         "otherwise.", call. = FALSE)
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
#' Flattens a [line_search()] object into the plain list the compiled loop
#' reads, so the C++ side needs no knowledge of S7. Every search produces the
#' same seven fields whichever subclass it is, and the compiled code branches
#' on `type` alone.
#'
#' @details
#' A field the given search has no use for is filled with a value the
#' compiled side ignores: `c2 = 0.9` for the two backtracking searches,
#' `shrink = 0.5` for Wolfe, `memory = 0` for both of the monotone ones.
#' Keeping the shape fixed is what lets one struct read all three.
#'
#' There are only **two** types. [nonmonotone()] describes itself as
#' `type = "armijo"` with `memory` above zero, the two differing in the
#' reference value alone, so the compiled loop needs no third branch. That is
#' the same fact the `memory = 0` identity records from the other side.
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
#' Armijo backtracking whose reference value is the worst of the last
#' `memory + 1` objective values instead of the current one, so a step may
#' make things worse now to be better placed later. Built by [nonmonotone()].
#'
#' @details
#' Beyond the `label` every line search carries, a `NonmonotoneSearch` holds
#' `c1`, `shrink`, `memory`, `max_step` and `resolution`. `memory` is the one
#' property no other line search has, and at `memory = 0` the object behaves
#' exactly as an [ArmijoSearch] with the same constants.
#'
#' @param c1 The sufficient-decrease constant.
#' @param shrink The factor the step is multiplied by on each backtrack.
#' @param memory How many earlier values the reference looks back over.
#' @param max_step The most backtracks allowed. A **count**.
#' @param resolution What the objective can tell apart, a number or a
#'   function of no arguments.
#'
#' @return An S7 object of class `NonmonotoneSearch` inheriting from
#'   [line_search()], with the five properties above beside `label`.
#'
#' @seealso [nonmonotone()] for the constructor, [bb()] for the method that
#'   needs it.
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
#' Armijo backtracking that compares against the worst of the last few
#' objective values rather than against the current one, so a step is allowed to
#' make things worse now in order to be better placed later.
#'
#' @param c1 Sufficient-decrease constant. Defaults to `1e-4`.
#' @param shrink Factor applied to the step on each backtrack. Defaults to
#'   `0.5`.
#' @param memory How many earlier values to look back over. Defaults to
#'   `10`; `0` makes this ordinary [armijo()].
#' @param max_step Most backtracks before the search gives up. Defaults to
#'   `30`.
#' @param resolution The smallest difference in the objective that means
#'   anything, in the objective's own units, or a function of no arguments
#'   returning it where it moves as the run goes. Defaults to `0`, which
#'   does not ask the question; see [armijo()] for what it is for and
#'   why it is asked at the full step rather than during the search.
#'
#' @details
#' The condition is Grippo, Lampariello and Lucidi's:
#' \deqn{f(x_k + s d_k) \le \max_{0 \le j \le m} f(x_{k-j}) + c_1 s\, g_k^\top d_k,}
#' which is Armijo's with the reference replaced by the largest of the last
#' \eqn{m+1} values. Every step it accepts improves on the worst of recent
#' memory; none is required to improve on the present.
#'
#' # Purpose
#'
#' Some methods are efficient *because* of steps that make the objective
#' worse. [bb()] is the clear case: its step length is a curvature estimate
#' taken from the last secant pair, and following that estimate faithfully
#' means occasionally going somewhere higher in order to be aligned with the
#' curvature when it matters. An Armijo condition forbids exactly those steps
#' and backtracks until it finds a shorter one, which is safe and is also most
#' of what the method was for.
#'
#' Measured on Rosenbrock, [bb()] takes 58 iterations and 67 objective
#' evaluations under `nonmonotone()` and 72 iterations and 154 evaluations
#' under [armijo()].
#'
#' The cost is that the guarantee weakens. A monotone method cannot cycle,
#' the objective being a decreasing sequence bounded below; a nonmonotone one
#' needs the finite memory to play that role, and the convergence result is
#' correspondingly more delicate. Use it where a method asks for it.
#'
#' # memory = 0
#'
#' At `memory = 0` the reference is the current value and the condition is
#' Armijo's exactly. Measured, `nonmonotone(memory = 0)` and [armijo()] give
#' the identical run, 72 iterations and 154 evaluations on the same problem,
#' so a comparison between the two settings is a comparison of the memory
#' alone.
#'
#' # No nonmonotone Wolfe
#'
#' The curvature condition is a statement about the gradient at the trial
#' point and says nothing about which value the decrease is measured against,
#' so a nonmonotone Wolfe search would be a fourth object rather than an
#' option on this one. There is not one here.
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
#' # What it buys the method it was added for.
#' f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
#' gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
#'                     200 * (p[2] - p[1]^2))
#' evals <- function(ls) minimize(bb(line_search = ls), f, c(-1.2, 1),
#'                                gr = gr)@counts[["f"]]
#' c(nonmonotone = evals(nonmonotone()), armijo = evals(armijo()))
#'
#' # And that memory = 0 is armijo(), which is what makes the comparison one
#' # of the memory alone.
#' identical(evals(nonmonotone(memory = 0)), evals(armijo()))
#'
#' @seealso [armijo()] for the monotone version, [bb()] for the method that
#'   needs this one, [line_search()] for the shared contract.
#' @export
nonmonotone <- function(c1 = 1e-4, shrink = 0.5, memory = 10, max_step = 30,
                        resolution = 0) {
  check_unit(c1, "c1")
  check_unit(shrink, "shrink")
  check_count(max_step, "max_step")
  check_resolution(resolution)
  if (length(memory) != 1L || !is.numeric(memory) || is.na(memory) ||
      memory < 0 || memory != round(memory)) {
    stop("'memory' must be a single non-negative whole number.", call. = FALSE)
  }
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
#' Checks that the value is a single number **strictly** inside
#' \eqn{(0, 1)}. Both endpoints are refused: `c1 = 0` asks for no decrease at
#' all and `c1 = 1` asks for the whole decrease the linear model predicts,
#' which a curved objective cannot supply.
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
#' Checks that the value is a single positive integer. Used for `max_step`,
#' which counts trials and not lengths, so a fractional value is a mistake
#' rather than a request.
#'
#' @param v The value.
#' @param nm Its name, for the message.
#'
#' @return Invisibly `TRUE`. Raises an error naming `nm` otherwise.
#'
#' @keywords internal
check_count <- function(v, nm) {
  if (length(v) != 1L || !is.numeric(v) || is.na(v) || v < 1) {
    stop("'", nm, "' must be a single positive number.", call. = FALSE)
  }
  invisible(TRUE)
}

# Fetched rather than captured, for the reason recorded in linkfunctions7:
# comparing S7 classes by identity breaks when the class object is re-created.
#' The line_search Class Object
#'
#' @description
#' Fetched rather than captured, so that a check cannot be fooled by the class
#' being re-created.
#'
#' @details
#' Comparing S7 classes by identity is object identity, so it is `FALSE`
#' for a class rebuilt from the same definition -- which is what happens under
#' any loader that re-evaluates the code rather than loading it, \pkg{covr}
#' among them. The same defect in \pkg{linkfunctions7} silently turned every
#' numerical fallback into a chain of first differences, and only the coverage
#' job noticed.
#'
#' @return The [line_search()] class object.
#'
#' @keywords internal
line_search_class <- function() line_search
