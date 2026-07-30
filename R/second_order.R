#' @include optimizer_class.R
#' @include generics.R
#' @include line_search.R
NULL

# Newton, BFGS and L-BFGS. Each is a constructor and a method that names its
# direction; the loop, the line search, the stopping rule and the reporting are
# all shared, which is what the frame was built for.


#' @title S7 Class for Newton's Method
#' @description The class \code{\link{newton}} instantiates.
#' @param hessian_mod How an indefinite Hessian is repaired.
#' @param floor The smallest eigenvalue the repaired Hessian may have.
#' @param step,line_search As in \code{\link{newton}}.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{newton}}
#' @name Newton-class
#' @aliases Newton
#' @keywords internal
Newton <- S7::new_class("Newton", parent = optimizer,
  properties = list(step = S7::class_numeric, line_search = S7::class_any,
                    hessian_mod = S7::class_character,
                    floor = S7::class_numeric))


#' @title Newton's Method with a Modified Hessian
#'
#' @description
#' Solves \eqn{H d = -g} for the direction, repairing \eqn{H} when it is not
#' positive definite, and then line searches along it.
#'
#' @param criterion The stopping rule; see \code{\link{crit_any}}.
#' @param hessian_mod How to repair an indefinite Hessian: \code{"eigen"}
#'   (default) or \code{"ridge"}. See Details.
#' @param floor The smallest eigenvalue the repaired Hessian is allowed.
#'   Defaults to \code{1e-8}.
#' @param step Initial step length offered to the line search. Defaults to
#'   \code{1}, which is the natural Newton step and is accepted unchanged near
#'   the solution.
#' @param line_search See \code{\link{armijo}} and \code{\link{wolfe}}.
#' @param maxit,max_eval,verbose,refresh,keep_trace As in
#'   \code{\link{gradient_descent}}.
#'
#' @details
#' Newton's method converges quadratically near a minimum and is unreliable away
#' from one, and the reason is worth stating plainly: \eqn{H^{-1}g} is only a
#' descent direction when \eqn{H} is positive definite. Where the objective
#' curves downwards the unmodified step points towards a saddle or a maximum, and
#' no line search can rescue it — every step along an ascent direction increases
#' the objective. This is not an edge case; it is the ordinary situation far from
#' the solution.
#'
#' Both repairs begin with an attempted Cholesky factorisation, which when it
#' succeeds is simultaneously the test for positive definiteness and the solve.
#' When it fails:
#'
#' \describe{
#'   \item{\code{"eigen"}}{decompose \eqn{H} and raise every eigenvalue below
#'     \code{floor} to it. The direction is then the Newton one in the subspace
#'     where the curvature is trustworthy, and gradient-like in the rest. Costs a
#'     symmetric eigendecomposition and gives the best-conditioned repair.}
#'   \item{\code{"ridge"}}{add \eqn{\tau I} with \eqn{\tau} doubling until the
#'     factorisation succeeds. This is Levenberg's idea: it interpolates between
#'     the Newton step at \eqn{\tau = 0} and a scaled steepest-descent step for
#'     large \eqn{\tau}. Cheaper, and blunter.}
#' }
#'
#' Which repair fired is recorded in the trace, so a run that spent its time
#' repairing rather than converging says so.
#'
#' \strong{On the Hessian itself.} If \code{he} is not supplied to
#' \code{\link{minimize}}, the Hessian is obtained by differencing the gradient.
#' That is one numerical differentiation when the gradient is analytic and
#' acceptable; when the gradient is \emph{also} differenced it is two composed,
#' which is the one place in the package where that happens and where the result
#' is correspondingly poor. An objective with neither derivative is better served
#' by \code{\link{bfgs}}, which never needs a Hessian at all.
#'
#' @return An S7 object of class \code{Newton}.
#'
#' @examples
#' rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' rosen_gr <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                           200 * (p[2] - p[1]^2))
#' rosen_he <- function(p) matrix(
#'   c(2 - 400 * (p[2] - 3 * p[1]^2), -400 * p[1],
#'     -400 * p[1], 200), 2, 2)
#'
#' minimize(newton(), rosen, c(-1.2, 1), gr = rosen_gr, he = rosen_he)
#'
#' @seealso \code{\link{bfgs}}, \code{\link{lbfgs}}, \code{\link{minimize}}
#' @export
newton <- function(criterion = crit_any(crit_grad(1e-8), crit_rel_obj(1e-12)),
                   hessian_mod = c("eigen", "ridge"), floor = 1e-8,
                   step = 1, line_search = armijo(),
                   maxit = 200, max_eval = 10000,
                   verbose = FALSE, refresh = 10, keep_trace = FALSE) {
  hessian_mod <- match.arg(hessian_mod)
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  if (length(floor) != 1L || !is.numeric(floor) || is.na(floor) || floor <= 0) {
    stop("'floor' must be a single positive number.", call. = FALSE)
  }
  Newton(name = "Newton", criterion = criterion, maxit = maxit,
         max_eval = max_eval, verbose = verbose, refresh = refresh,
         keep_trace = keep_trace, step = step, line_search = line_search,
         hessian_mod = hessian_mod, floor = floor)
}


#' @title S7 Class for BFGS
#' @description The class \code{\link{bfgs}} instantiates.
#' @param curv_tol The curvature threshold below which the update is skipped.
#' @param max_skip Consecutive skips before the approximation is reset.
#' @param step,line_search As in \code{\link{bfgs}}.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{bfgs}}
#' @name Bfgs-class
#' @aliases Bfgs
#' @keywords internal
Bfgs <- S7::new_class("Bfgs", parent = optimizer,
  properties = list(step = S7::class_numeric, line_search = S7::class_any,
                    curv_tol = S7::class_numeric,
                    max_skip = S7::class_numeric))


#' @title BFGS
#'
#' @description
#' Builds an approximation to the inverse Hessian from successive gradients, so
#' the direction is a matrix-vector product and no second derivatives are ever
#' required.
#'
#' @param criterion The stopping rule.
#' @param curv_tol The update is skipped when
#'   \eqn{s^\top y \le \texttt{curv\_tol}\,\lVert s\rVert\lVert y\rVert}.
#'   Defaults to \code{1e-10}.
#' @param max_skip Consecutive skipped updates before the approximation is reset
#'   to the identity. Defaults to 5.
#' @param step,line_search,maxit,max_eval,verbose,refresh,keep_trace As in
#'   \code{\link{newton}}. The default line search is \code{\link{wolfe}}; see
#'   Details.
#'
#' @details
#' The default line search is the strong Wolfe one, and that is a requirement
#' rather than a preference. BFGS updates its approximation from the secant pair
#' \eqn{(s, y)} with \eqn{s = x_{new} - x_{old}} and
#' \eqn{y = g_{new} - g_{old}}, and the update is only meaningful when
#' \eqn{s^\top y > 0}. Armijo backtracking can accept a step so short that the
#' gradient has barely moved, giving a pair with no curvature information in it;
#' the Wolfe curvature condition is exactly the guarantee that this does not
#' happen. Run on Armijo alone, BFGS loses the property that makes it BFGS.
#'
#' When the curvature condition fails anyway the update is \strong{skipped}
#' rather than applied. A small \eqn{s^\top y} makes \eqn{\rho = 1/s^\top y}
#' enormous and one bad step destroys the accumulated approximation; a stale but
#' sound matrix is better than a fresh but corrupted one. After
#' \code{max_skip} consecutive skips there is no curvature information left worth
#' keeping, and the matrix is reset to the identity — the method restarts as
#' steepest descent and rebuilds.
#'
#' The first accepted pair also rescales the identity by
#' \eqn{s^\top y / y^\top y}. Without it the first quasi-Newton step is taken
#' with a unit Hessian, which on a badly scaled problem has entirely the wrong
#' magnitude and wastes a line search discovering so.
#'
#' @return An S7 object of class \code{Bfgs}.
#'
#' @examples
#' rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' minimize(bfgs(), rosen, c(-1.2, 1))
#'
#' @seealso \code{\link{lbfgs}} for many parameters, \code{\link{newton}} when a
#'   Hessian is available.
#' @export
bfgs <- function(criterion = crit_any(crit_grad(1e-8), crit_rel_obj(1e-12)),
                 curv_tol = 1e-10, max_skip = 5,
                 step = 1, line_search = wolfe(),
                 maxit = 500, max_eval = 20000,
                 verbose = FALSE, refresh = 10, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  Bfgs(name = "BFGS", criterion = criterion, maxit = maxit,
       max_eval = max_eval, verbose = verbose, refresh = refresh,
       keep_trace = keep_trace, step = step, line_search = line_search,
       curv_tol = curv_tol, max_skip = max_skip)
}


#' @title S7 Class for Limited-Memory BFGS
#' @description The class \code{\link{lbfgs}} instantiates.
#' @param memory How many secant pairs to keep.
#' @param curv_tol The curvature threshold below which a pair is not stored.
#' @param step,line_search As in \code{\link{lbfgs}}.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{lbfgs}}
#' @name Lbfgs-class
#' @aliases Lbfgs
#' @keywords internal
Lbfgs <- S7::new_class("Lbfgs", parent = optimizer,
  properties = list(step = S7::class_numeric, line_search = S7::class_any,
                    memory = S7::class_numeric,
                    curv_tol = S7::class_numeric))


#' @title Limited-Memory BFGS
#'
#' @description
#' BFGS without ever forming the matrix: only the last \code{memory} secant
#' pairs are kept, and the direction comes from the two-loop recursion.
#'
#' @param criterion The stopping rule.
#' @param memory How many secant pairs to keep. Defaults to 10.
#' @param curv_tol A pair is discarded when
#'   \eqn{s^\top y \le \texttt{curv\_tol}\,\lVert s\rVert\lVert y\rVert}.
#'   Defaults to \code{1e-10}.
#' @param step,line_search,maxit,max_eval,verbose,refresh,keep_trace As in
#'   \code{\link{bfgs}}.
#'
#' @details
#' The cost is \eqn{O(mp)} in both time and memory rather than \eqn{O(p^2)}. For
#' a handful of parameters that is no gain whatever and \code{\link{bfgs}} will
#' converge in fewer iterations, since it carries the whole approximation; the
#' crossover is somewhere in the hundreds of parameters, and beyond a few
#' thousand the full matrix is simply not storable.
#'
#' The recursion is scaled at each iteration by the most recent pair's
#' \eqn{s^\top y / y^\top y}. That single number is doing the work the full
#' matrix would otherwise do, and it is the reason the method converges at all
#' without one.
#'
#' A larger \code{memory} is not uniformly better: old pairs describe curvature
#' at points the iterate has left, and on a strongly nonlinear objective they can
#' be worse than no information. Ten is the conventional choice for good reason.
#'
#' @return An S7 object of class \code{Lbfgs}.
#'
#' @examples
#' rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' minimize(lbfgs(memory = 5), rosen, c(-1.2, 1))
#'
#' @seealso \code{\link{bfgs}}
#' @export
lbfgs <- function(criterion = crit_any(crit_grad(1e-8), crit_rel_obj(1e-12)),
                  memory = 10, curv_tol = 1e-10,
                  step = 1, line_search = wolfe(),
                  maxit = 500, max_eval = 20000,
                  verbose = FALSE, refresh = 10, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  check_count(memory, "memory")
  Lbfgs(name = "L-BFGS", criterion = criterion, maxit = maxit,
        max_eval = max_eval, verbose = verbose, refresh = refresh,
        keep_trace = keep_trace, step = step, line_search = line_search,
        memory = memory, curv_tol = curv_tol)
}


# --- the methods ------------------------------------------------------------
#
# Each names its direction and hands everything else to run_descent(); that they
# are three lines apiece is the point of the frame.

#' @title Minimise by Newton's Method
#' @name minimize.Newton
#' @param optimizer A \code{Newton} object.
#' @param fn,par,gr,he,bounds,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Newton) <-
  function(optimizer, fn, par, gr = NULL, he = NULL, bounds = NULL, ...) {
    run_descent(optimizer, fn, par, gr, he, bounds,
                list(type = "newton", hessian_mod = optimizer@hessian_mod,
                     floor = optimizer@floor))
  }

#' @title Minimise by BFGS
#' @name minimize.Bfgs
#' @param optimizer A \code{Bfgs} object.
#' @param fn,par,gr,he,bounds,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Bfgs) <-
  function(optimizer, fn, par, gr = NULL, he = NULL, bounds = NULL, ...) {
    run_descent(optimizer, fn, par, gr, he, bounds,
                list(type = "bfgs", curv_tol = optimizer@curv_tol,
                     max_skip = as.integer(optimizer@max_skip)))
  }

#' @title Minimise by Limited-Memory BFGS
#' @name minimize.Lbfgs
#' @param optimizer An \code{Lbfgs} object.
#' @param fn,par,gr,he,bounds,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Lbfgs) <-
  function(optimizer, fn, par, gr = NULL, he = NULL, bounds = NULL, ...) {
    run_descent(optimizer, fn, par, gr, he, bounds,
                list(type = "lbfgs", memory = as.integer(optimizer@memory),
                     curv_tol = optimizer@curv_tol))
  }
