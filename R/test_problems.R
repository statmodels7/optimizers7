#' @include criterion.R
NULL

#' @title The Standard Test Problems
#'
#' @description
#' Returns the eight functions optimization papers are argued over: a
#' quadratic, a curved valley, a few awkward polynomials, two with many
#' minima, and one with a kink. Each carries its analytic gradient, a
#' customary starting point, its known minimizer and two flags describing
#' the difficulty it poses. They are exported for testing any optimizer, not
#' only this package's; [check_optimizer()] runs the whole battery.
#'
#' @details
#' # The eight
#'
#' \tabular{llll}{
#'   **name** \tab **p** \tab **start** \tab **what it tests** \cr
#'   `sphere` \tab 3 \tab (1.3, -0.7, 0.8) \tab nothing; a run that fails here is broken \cr
#'   `rosenbrock` \tab 2 \tab (-1.2, 1) \tab a curved valley, Hessian condition 2510 at the solution \cr
#'   `booth` \tab 2 \tab (0, 0) \tab a well-conditioned quadratic in disguise \cr
#'   `beale` \tab 2 \tab (1, 1) \tab a narrow valley, Hessian condition 163 \cr
#'   `powell` \tab 4 \tab (3, -1, 0, 1) \tab a **singular** Hessian at the solution \cr
#'   `himmelblau` \tab 2 \tab (0, 0) \tab four minima, all of value zero \cr
#'   `rastrigin` \tab 2 \tab (0.4, -0.4) \tab 121 local minima on the usual box \cr
#'   `abs_sum` \tab 3 \tab (0, 0, 0) \tab a kink at the solution itself \cr
#' }
#'
#' The starting points are the customary ones, which for `rosenbrock` and
#' `powell` are the deliberately unhelpful points those functions were
#' designed around. The two hardest are hard for opposite reasons.
#' `powell`'s Hessian at the origin has eigenvalues 202, 20, 0 and 0, so a
#' second-order method has no curvature to read in two of its four
#' directions and converges linearly where it usually converges
#' quadratically. `rastrigin` has curvature everywhere and 121 places to
#' stop, so what decides the answer is where the run began.
#'
#' # The two flags
#'
#' `multimodal` marks a surface with more than one local minimum. A local
#' method that reaches a different one there is behaving correctly, so
#' scoring it against `solution` would report a failure that is not one.
#' [check_optimizer()] therefore scores every problem on the value reached
#' and labels these two in its `note` column.
#'
#' `smooth` is `FALSE` for `abs_sum` alone. Its gradient at the solution is
#' exactly zero, because `sign(0)` is zero, and at every neighboring point
#' the subgradient has max-norm one. A method that tests \eqn{\nabla f} will
#' therefore arrive at the answer and be unable to certify it; the
#' derivative-free and non-smooth methods report a measure of their own
#' instead, tested by [crit_stationary()].
#'
#' # Every solution has value zero, and that hides an effect
#'
#' A line search accepts a step only when the objective decreases by a
#' definite amount, and near a minimum that decrease is about
#' \eqn{\lVert g \rVert^{2} / (2\lambda)} for a curvature \eqn{\lambda}.
#' Once it falls below the rounding of the objective itself, about
#' \eqn{\varepsilon \lvert f \rvert}, no step in any direction can be
#' verified and the search refuses all of them. The smallest gradient a run
#' can reach is therefore of order
#'
#' \deqn{\lVert g \rVert_{\text{floor}} \approx
#'   \sqrt{2 \lambda \varepsilon \lvert f(x^{\ast}) \rvert},}
#'
#' which grows with the value at the solution and is exactly zero for every
#' problem here. A log-likelihood is the opposite case, of order one at its
#' optimum, so an optimizer that reaches `1e-15` on this battery may stop at
#' `1e-8` there. The defaults of [crit_grad()] allow for that.
#'
#' Adding a constant to any of these objectives moves neither the minimizer
#' nor the gradient, and reproduces the effect on demand. Conjugate
#' gradients on `rosenbrock`, asked for `crit_grad(1e-14)` with a budget of
#' 20000 iterations, reaches a max-norm gradient of `4.4e-15` at
#' \eqn{f^{\ast} = 0} and stops at `4.4e-08`, `2.8e-06` and `6.5e-05` once
#' the constants \eqn{1}, \eqn{10^{3}} and \eqn{10^{6}} are added. Only the
#' first run reports `converged = TRUE`; in the other three the tolerance
#' asked for was never reachable.
#'
#' @param which A character vector of problem names, or `NULL`. `NULL`, the
#'   default, returns all eight in the order of the table above; a character
#'   vector returns just those, in the order given, so
#'   `test_problems(c("beale", "sphere"))` comes back beale first. A name
#'   that is not one of the eight raises an error listing the eight.
#'
#' @return A named list, one element per problem, each itself a list of eight
#'   components:
#'   \describe{
#'     \item{`name`}{character, the same as the element's name.}
#'     \item{`fn`}{the objective, a function of a numeric vector of length
#'       `p` returning one number.}
#'     \item{`gr`}{its analytic gradient, returning a vector of length `p`.}
#'     \item{`par`}{numeric of length `p`, the customary starting point.}
#'     \item{`solution`}{numeric of length `p`, the minimizer. For
#'       `himmelblau` this is the one of the four nearest the start.}
#'     \item{`value`}{numeric, the minimum. Zero for all eight.}
#'     \item{`multimodal`}{logical, `TRUE` for `himmelblau` and
#'       `rastrigin`.}
#'     \item{`smooth`}{logical, `FALSE` for `abs_sum` alone.}
#'   }
#'
#' @examples
#' names(test_problems())
#'
#' # Every problem's stated solution really is a stationary point of value zero.
#' P <- test_problems()
#' stopifnot(all(vapply(P, function(q) q$fn(q$solution), 0) == 0))
#'
#' # BFGS solves the curved valley from the unhelpful start.
#' p <- test_problems("rosenbrock")[[1]]
#' fit <- minimize(bfgs(), p$fn, p$par, gr = p$gr)
#' all.equal(fit@par, p$solution, tolerance = 1e-6)
#'
#' # Powell's difficulty is a singular Hessian at the solution: the quartic
#' # terms contribute nothing to the curvature at the origin, so two of the
#' # four eigenvalues are exactly zero and a second-order method has no
#' # curvature to read in those directions.
#' H <- matrix(c(2, 20, 0, 0,
#'               20, 200, 0, 0,
#'               0, 0, 10, -10,
#'               0, 0, -10, 10), 4, 4)
#' eigen(H, only.values = TRUE)$values
#'
#' # abs_sum has a zero gradient at the solution and a subgradient of max-norm
#' # one at every neighboring point, so a gradient rule cannot certify it.
#' k <- test_problems("abs_sum")[[1]]
#' k$gr(k$solution)
#' k$gr(k$solution + 1e-9)
#'
#' # A name that is not one of the eight is refused, and the message lists them.
#' try(test_problems("banana"))
#'
#' @seealso [check_optimizer()], which runs the battery and reports what each
#'   method reached, and [crit_stationary()] for the rule the non-smooth
#'   problem needs.
#' @export
test_problems <- function(which = NULL) {
  P <- list(

    sphere = list(
      name = "sphere",
      fn = function(p) sum(p^2),
      gr = function(p) 2 * p,
      par = c(1.3, -0.7, 0.8),
      solution = c(0, 0, 0), value = 0,
      multimodal = FALSE, smooth = TRUE),

    rosenbrock = list(
      name = "rosenbrock",
      fn = function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2,
      gr = function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
                         200 * (p[2] - p[1]^2)),
      par = c(-1.2, 1),
      solution = c(1, 1), value = 0,
      multimodal = FALSE, smooth = TRUE),

    booth = list(
      name = "booth",
      fn = function(p) (p[1] + 2 * p[2] - 7)^2 + (2 * p[1] + p[2] - 5)^2,
      gr = function(p) {
        a <- p[1] + 2 * p[2] - 7; b <- 2 * p[1] + p[2] - 5
        c(2 * a + 4 * b, 4 * a + 2 * b)
      },
      par = c(0, 0),
      solution = c(1, 3), value = 0,
      multimodal = FALSE, smooth = TRUE),

    beale = list(
      name = "beale",
      fn = function(p) {
        x <- p[1]; y <- p[2]
        (1.5 - x + x * y)^2 + (2.25 - x + x * y^2)^2 + (2.625 - x + x * y^3)^2
      },
      gr = function(p) {
        x <- p[1]; y <- p[2]
        t1 <- 1.5 - x + x * y; t2 <- 2.25 - x + x * y^2
        t3 <- 2.625 - x + x * y^3
        c(2 * t1 * (y - 1) + 2 * t2 * (y^2 - 1) + 2 * t3 * (y^3 - 1),
          2 * t1 * x + 4 * t2 * x * y + 6 * t3 * x * y^2)
      },
      par = c(1, 1),
      solution = c(3, 0.5), value = 0,
      multimodal = FALSE, smooth = TRUE),

    powell = list(
      name = "powell",
      fn = function(p) {
        (p[1] + 10 * p[2])^2 + 5 * (p[3] - p[4])^2 +
          (p[2] - 2 * p[3])^4 + 10 * (p[1] - p[4])^4
      },
      gr = function(p) {
        a <- p[1] + 10 * p[2]; b <- p[3] - p[4]
        cc <- p[2] - 2 * p[3]; d <- p[1] - p[4]
        c(2 * a + 40 * d^3,
          20 * a + 4 * cc^3,
          10 * b - 8 * cc^3,
          -10 * b - 40 * d^3)
      },
      par = c(3, -1, 0, 1),
      solution = c(0, 0, 0, 0), value = 0,
      multimodal = FALSE, smooth = TRUE),

    himmelblau = list(
      name = "himmelblau",
      fn = function(p) (p[1]^2 + p[2] - 11)^2 + (p[1] + p[2]^2 - 7)^2,
      gr = function(p) {
        a <- p[1]^2 + p[2] - 11; b <- p[1] + p[2]^2 - 7
        c(4 * p[1] * a + 2 * b, 2 * a + 4 * p[2] * b)
      },
      par = c(0, 0),
      # Four minima, all with value zero; (3, 2) is the one nearest the start.
      solution = c(3, 2), value = 0,
      multimodal = TRUE, smooth = TRUE),

    rastrigin = list(
      name = "rastrigin",
      fn = function(p) 10 * length(p) + sum(p^2 - 10 * cos(2 * pi * p)),
      gr = function(p) 2 * p + 20 * pi * sin(2 * pi * p),
      par = c(0.4, -0.4),
      solution = c(0, 0), value = 0,
      multimodal = TRUE, smooth = TRUE),

    abs_sum = list(
      name = "abs_sum",
      fn = function(p) sum(abs(p - c(1, -2, 0.5))),
      gr = function(p) sign(p - c(1, -2, 0.5)),
      par = c(0, 0, 0),
      solution = c(1, -2, 0.5), value = 0,
      multimodal = FALSE, smooth = FALSE)
  )

  if (is.null(which)) return(P)
  missing <- setdiff(which, names(P))
  if (length(missing)) {
    stop("No such test problem: ", paste(missing, collapse = ", "),
         ". Available: ", paste(names(P), collapse = ", "), call. = FALSE)
  }
  P[which]
}
