#' @include criterion.R
NULL

#' @title The Standard Test Problems
#'
#' @description
#' The functions optimisation papers are argued over: a quadratic, a curved
#' valley, a few awkward polynomials, two with many minima, and one with a kink.
#' Each carries its analytic gradient and its known answer.
#'
#' @param which An optional character vector naming a subset.
#'
#' @details
#' They are exported because they are useful to anyone writing an optimiser, not
#' only to this package's own tests: \code{\link{check_optimizer}} runs them, and
#' so can you.
#'
#' Each element is a list with \code{name}, \code{fn}, \code{gr}, \code{par}
#' (a starting point), \code{solution}, \code{value}, and two flags.
#' \code{multimodal} marks a surface with more than one local minimum, where a
#' local method reaching a different one is behaving correctly and not failing;
#' \code{smooth} is \code{FALSE} for the one whose derivative does not exist
#' everywhere, where a method that assumes it does will arrive at the answer and
#' then be unable to certify it.
#'
#' The starting points are the ones customarily used, which for
#' \code{rosenbrock} and \code{powell} means the deliberately unhelpful ones the
#' functions were designed around.
#'
#' @return A named list of problems.
#'
#' @examples
#' names(test_problems())
#'
#' p <- test_problems("rosenbrock")[[1]]
#' minimize(bfgs(), p$fn, p$par, gr = p$gr)@par
#' p$solution
#'
#' @seealso \code{\link{check_optimizer}}
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
