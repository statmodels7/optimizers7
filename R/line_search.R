#' @include criterion.R
NULL

#' @title S7 Class for Line Searches
#'
#' @description
#' How far to step along a direction, as an object.
#'
#' @details
#' A descent method produces a \emph{direction}; how far to travel along it is a
#' separate question with its own theory. Separating the two is what lets Newton,
#' BFGS and L-BFGS share one carefully written answer instead of each carrying
#' its own, and what lets a user change the answer without touching the method.
#'
#' The class is abstract; use \code{\link{armijo}} or \code{\link{wolfe}}.
#'
#' @param label A short character label, used when printing.
#'
#' @return An S7 object of class \code{line_search}.
#'
#' @examples
#' # Abstract: use one of the constructors.
#' try(line_search(label = "mine"))
#'
#' armijo()
#' wolfe()
#'
#' # A method takes whichever it is given, and the choice shows in the result.
#' f <- function(p) sum((p - c(1, 2))^2)
#' g <- function(p) 2 * (p - c(1, 2))
#' minimize(bfgs(line_search = armijo()), f, c(0, 0), gr = g)@counts
#' minimize(bfgs(line_search = wolfe()), f, c(0, 0), gr = g)@counts
#'
#' @seealso \code{\link{armijo}}, \code{\link{wolfe}}
#' @export
line_search <- S7::new_class(
  "line_search",
  properties = list(label = S7::class_character),
  abstract = TRUE
)


#' @title S7 Class for Armijo Backtracking
#' @description The class \code{\link{armijo}} instantiates.
#' @param c1 The sufficient-decrease constant.
#' @param shrink The factor the step is multiplied by on each backtrack.
#' @param max_step The most backtracks allowed.
#' @return An S7 object inheriting from \code{\link{line_search}}.
#' @seealso \code{\link{armijo}}
#' @keywords internal
ArmijoSearch <- S7::new_class("ArmijoSearch", parent = line_search,
  properties = list(c1 = S7::class_numeric, shrink = S7::class_numeric,
                    max_step = S7::class_numeric))


#' @title Backtracking Line Search with the Armijo Condition
#'
#' @description
#' Shrink the step until the objective decreases by enough:
#' \deqn{f(x + s d) \leq f(x) + c_1 s\, g^\top d .}
#'
#' @param c1 Sufficient-decrease constant, in \eqn{(0, 1)}. Defaults to
#'   \code{1e-4}, the conventional value: it demands a decrease, but only a tiny
#'   fraction of what the linear model predicts, so it almost never rejects a
#'   sensible step.
#' @param shrink Factor applied on each backtrack, in \eqn{(0, 1)}. Defaults to
#'   \code{0.5}.
#' @param max_step Maximum backtracks before the search gives up. Defaults to 30.
#'
#' @details
#' The \eqn{c_1 s\, g^\top d} term is the whole point, and dropping it — testing
#' merely \eqn{f_{new} \le f} — is a real defect rather than a simplification.
#' On a quadratic with unit step the gradient update reflects the iterate through
#' the minimum, leaving the objective \emph{exactly} unchanged; the weak test
#' accepts it, the iterate oscillates forever, and a stopping rule watching the
#' objective sees no change and reports convergence at a point that is not a
#' minimum.
#'
#' Cheap: it evaluates the objective at trial points and never the gradient. That
#' is enough for a method that only needs to make progress, and not enough for a
#' quasi-Newton method, which needs \code{\link{wolfe}}.
#'
#' @return A \code{\link{line_search}} object.
#'
#' @examples
#' armijo()
#' minimize(gd(line_search = armijo(shrink = 0.2)),
#'          function(p) sum((p - 1:2)^2), c(0, 0))
#'
#' @seealso \code{\link{wolfe}}
#' @export
armijo <- function(c1 = 1e-4, shrink = 0.5, max_step = 30) {
  check_unit(c1, "c1")
  check_unit(shrink, "shrink")
  check_count(max_step, "max_step")
  ArmijoSearch(label = paste0("Armijo backtracking (c1 = ", format(c1), ")"),
               c1 = c1, shrink = shrink, max_step = max_step)
}


#' @title S7 Class for the Strong Wolfe Line Search
#' @description The class \code{\link{wolfe}} instantiates.
#' @param c1 The sufficient-decrease constant.
#' @param c2 The curvature constant.
#' @param max_step The most trial steps allowed.
#' @return An S7 object inheriting from \code{\link{line_search}}.
#' @seealso \code{\link{wolfe}}
#' @keywords internal
WolfeSearch <- S7::new_class("WolfeSearch", parent = line_search,
  properties = list(c1 = S7::class_numeric, c2 = S7::class_numeric,
                    max_step = S7::class_numeric))


#' @title Line Search Satisfying the Strong Wolfe Conditions
#'
#' @description
#' Sufficient decrease and a curvature condition together:
#' \deqn{f(x + s d) \leq f(x) + c_1 s\, g^\top d, \qquad
#'       \lvert g(x + s d)^\top d \rvert \leq c_2 \lvert g^\top d \rvert .}
#'
#' @param c1 Sufficient-decrease constant. Defaults to \code{1e-4}.
#' @param c2 Curvature constant, with \eqn{c_1 < c_2 < 1}. Defaults to
#'   \code{0.9}, the usual choice for a quasi-Newton method; \code{0.1} is
#'   conventional for nonlinear conjugate gradients, which want a more exact
#'   line search.
#' @param max_step Maximum trial steps in each of the bracketing and zoom
#'   phases. Defaults to 30.
#'
#' @details
#' The curvature condition is what \code{\link{armijo}} cannot provide, and it is
#' not a refinement: a quasi-Newton method needs it to work at all. BFGS builds
#' its approximation from the secant pair \eqn{(s, y)} with
#' \eqn{y = g_{new} - g_{old}}, and a step so short that the gradient has barely
#' moved gives a pair carrying no curvature information — so the update must
#' either be skipped or it corrupts the matrix. Requiring the gradient along the
#' direction to have shrunk by a factor \eqn{c_2} is exactly the guarantee that
#' this does not happen.
#'
#' The implementation is the bracketing-and-zoom scheme, with bisection inside
#' the zoom rather than polynomial interpolation: a few more evaluations, and it
#' cannot be defeated by an awkwardly shaped interval.
#'
#' It costs gradient evaluations at trial points, which \code{\link{armijo}} does
#' not, so it is the more expensive choice per iteration and usually the cheaper
#' one per problem.
#'
#' @return A \code{\link{line_search}} object.
#'
#' @examples
#' wolfe()
#' minimize(gd(line_search = wolfe(), maxit = 200),
#'          function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2, c(-1.2, 1))
#'
#' @seealso \code{\link{armijo}}
#' @export
wolfe <- function(c1 = 1e-4, c2 = 0.9, max_step = 30) {
  check_unit(c1, "c1")
  check_unit(c2, "c2")
  if (c2 <= c1) {
    stop("'c2' must be greater than 'c1'; the conditions are unsatisfiable ",
         "otherwise.", call. = FALSE)
  }
  check_count(max_step, "max_step")
  WolfeSearch(label = paste0("strong Wolfe (c1 = ", format(c1),
                             ", c2 = ", format(c2), ")"),
              c1 = c1, c2 = c2, max_step = max_step)
}


#' Describe a Line Search to the C++ Side
#'
#' @description
#' Flattens a \code{\link{line_search}} object into the list the compiled code
#' reads.
#'
#' @param x A \code{\link{line_search}} object.
#'
#' @return A list with \code{type}, \code{c1}, \code{c2}, \code{shrink} and
#'   \code{max_step}. Fields a given search does not use are filled with values
#'   the C++ side ignores, so that the structure is the same shape whichever
#'   search it describes.
#'
#' @keywords internal
line_search_spec <- S7::new_generic("line_search_spec", "x",
                                    function(x) S7::S7_dispatch())

S7::method(line_search_spec, ArmijoSearch) <- function(x) {
  list(type = "armijo", c1 = x@c1, c2 = 0.9, shrink = x@shrink,
       max_step = as.integer(x@max_step), memory = 0L)
}

S7::method(line_search_spec, WolfeSearch) <- function(x) {
  list(type = "wolfe", c1 = x@c1, c2 = x@c2, shrink = 0.5,
       max_step = as.integer(x@max_step), memory = 0L)
}


#' @title S7 Class for the Nonmonotone Line Search
#' @description The class \code{\link{nonmonotone}} instantiates.
#' @param c1 The sufficient-decrease constant.
#' @param shrink The factor the step is multiplied by on each backtrack.
#' @param memory How many earlier values the reference looks back over.
#' @param max_step The most backtracks allowed.
#' @return An S7 object inheriting from \code{\link{line_search}}.
#' @seealso \code{\link{nonmonotone}}
#' @keywords internal
NonmonotoneSearch <- S7::new_class("NonmonotoneSearch", parent = line_search,
  properties = list(
    c1       = S7::class_numeric,
    shrink   = S7::class_numeric,
    memory   = S7::class_numeric,
    max_step = S7::class_numeric
  ))


#' @title Nonmonotone Backtracking
#'
#' @description
#' Armijo backtracking that compares against the worst of the last few
#' objective values rather than against the current one, so a step is allowed to
#' make things worse now in order to be better placed later.
#'
#' @param c1 Sufficient-decrease constant. Defaults to \code{1e-4}.
#' @param shrink Factor applied to the step on each backtrack. Defaults to
#'   \code{0.5}.
#' @param memory How many earlier values to look back over. Defaults to
#'   \code{10}; \code{0} makes this ordinary \code{\link{armijo}}.
#' @param max_step Most backtracks before the search gives up. Defaults to
#'   \code{30}.
#'
#' @details
#' The condition is Grippo, Lampariello and Lucidi's:
#' \deqn{f(x_k + s d_k) \le \max_{0 \le j \le m} f(x_{k-j}) + c_1 s\, g_k^	op d_k,}
#' which is Armijo's with the reference replaced by the largest of the last
#' \eqn{m+1} values. Every step it accepts improves on the worst of recent
#' memory; none is required to improve on the present.
#'
#' \subsection{What it is for}{
#' Some methods are efficient \emph{because} of steps that make the objective
#' worse. \code{\link{bb}} is the clear case: its step length is a curvature
#' estimate taken from the last secant pair, and following that estimate
#' faithfully means occasionally going somewhere higher in order to be aligned
#' with the curvature when it matters. An Armijo condition forbids exactly those
#' steps and backtracks until it finds a shorter one, which is safe and is also
#' most of what the method was for.
#'
#' The cost is that the guarantee weakens. A monotone method cannot cycle,
#' because the objective is a decreasing sequence bounded below; a nonmonotone
#' one needs the finite memory to play that role, and the convergence result is
#' correspondingly more delicate. Use it where a method asks for it, not as a
#' faster default.
#' }
#'
#' \subsection{Not available with Wolfe}{
#' The curvature condition is a statement about the gradient at the trial point
#' and has nothing to say about which value the decrease is measured against, so
#' a nonmonotone Wolfe search would be a different object rather than an option
#' on this one. There is not one here.
#' }
#'
#' @return A \code{\link{line_search}} object.
#'
#' @references
#' Grippo, L., Lampariello, F. and Lucidi, S. (1986). A nonmonotone line search
#' technique for Newton's method. \emph{SIAM Journal on Numerical Analysis}
#' \strong{23}, 707--716.
#'
#' Raydan, M. (1997). The Barzilai and Borwein gradient method for the large
#' scale unconstrained minimization problem. \emph{SIAM Journal on
#' Optimization} \strong{7}, 26--33.
#'
#' @examples
#' nonmonotone()
#' nonmonotone(memory = 5)
#'
#' # what it buys the method it was added for
#' f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
#' gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
#'                     200 * (p[2] - p[1]^2))
#' minimize(bb(line_search = armijo()), f, c(-1.2, 1), gr = gr)@iterations
#' minimize(bb(), f, c(-1.2, 1), gr = gr)@iterations
#'
#' @seealso \code{\link{armijo}}, \code{\link{bb}}
#' @export
nonmonotone <- function(c1 = 1e-4, shrink = 0.5, memory = 10, max_step = 30) {
  check_unit(c1, "c1")
  check_unit(shrink, "shrink")
  check_count(max_step, "max_step")
  if (length(memory) != 1L || !is.numeric(memory) || is.na(memory) ||
      memory < 0 || memory != round(memory)) {
    stop("'memory' must be a single non-negative whole number.", call. = FALSE)
  }
  NonmonotoneSearch(
    label = paste0("nonmonotone backtracking (memory = ", format(memory), ")"),
    c1 = c1, shrink = shrink, memory = memory, max_step = max_step)
}


#' @rdname line_search_spec
#' @name line_search_spec.NonmonotoneSearch
#' @keywords internal
S7::method(line_search_spec, NonmonotoneSearch) <- function(x) {
  list(type = "armijo", c1 = x@c1, c2 = 0.9, shrink = x@shrink,
       max_step = as.integer(x@max_step), memory = as.integer(x@memory))
}

#' @title Print Method for Line Searches
#' @name print.line_search
#' @param x A \code{\link{line_search}} object.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @examples
#' print(wolfe())
#' @keywords internal
S7::method(print, line_search) <- function(x, ...) {
  cat("<line_search> ", x@label, "\n", sep = "")
  invisible(x)
}


#' Validate a Constant in the Unit Interval
#'
#' @param v The value.
#' @param nm Its name, for the message.
#'
#' @return Invisibly \code{TRUE}; raises an error otherwise.
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
#' @param v The value.
#' @param nm Its name, for the message.
#'
#' @return Invisibly \code{TRUE}; raises an error otherwise.
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
#' Comparing S7 classes by identity is object identity, so it is \code{FALSE}
#' for a class rebuilt from the same definition -- which is what happens under
#' any loader that re-evaluates the code rather than loading it, \pkg{covr}
#' among them. The same trap in \pkg{linkfunctions7} silently turned every
#' numerical fallback into a chain of first differences, and only the coverage
#' job noticed.
#'
#' @return The \code{\link{line_search}} class object.
#'
#' @keywords internal
line_search_class <- function() line_search
