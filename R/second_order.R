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
#'   \code{\link{gd}}.
#'
#' @details
#' Newton's method converges quadratically near a minimum and is unreliable away
#' from one: \eqn{H^{-1}g} is only a
#' descent direction when \eqn{H} is positive definite. Where the objective
#' curves downwards the unmodified step points towards a saddle or a maximum, and
#' no line search can rescue it — every step along an ascent direction increases
#' the objective. This is not an edge case; it is the ordinary situation far from
#' the solution.
#'
#' Both repairs begin with an attempted Cholesky factorization, which when it
#' succeeds is simultaneously the test for positive definiteness and the solve.
#' When it fails:
#'
#' \describe{
#'   \item{\code{"eigen"}}{decompose \eqn{H} and raise every eigenvalue below
#'     \code{floor} to it. The direction is then the Newton one in the subspace
#'     where the curvature is trustworthy, and gradient-like in the rest. Costs a
#'     symmetric eigendecomposition and gives the best-conditioned repair.}
#'   \item{\code{"ridge"}}{add \eqn{\tau I} with \eqn{\tau} doubling until the
#'     factorization succeeds. This is Levenberg's idea: it interpolates between
#'     the Newton step at \eqn{\tau = 0} and a scaled steepest-descent step for
#'     large \eqn{\tau}. Cheaper, and blunter.}
#' }
#'
#' Which repair fired is recorded in the trace, so the trace shows a run
#' that spent its time repairing rather than converging.
#'
#' \strong{A repaired step is capped, and an unrepaired one is not.} Where the
#' factorization succeeds, \eqn{H^{-1}g} is the Newton step and \eqn{step = 1}
#' is its own unit, so it is taken as it is. Where it fails, the component of
#' the direction in the floored subspace is \eqn{g_i/\lambda_{\mathrm{floor}}},
#' whose size is set by \code{floor} rather than by any curvature of the
#' objective: the direction is therefore scaled to
#' \eqn{\min(1, 1/\lVert d\rVert_\infty)}, a displacement of order one in the
#' parameters. This is the rule \code{\link{bb}} applies when its secant pair
#' reports no curvature and the one \code{\link{bfgs}} and \code{\link{lbfgs}}
#' apply to their first direction, and it ONLY EVER SHORTENS, so a run whose
#' repaired steps were already of order one is unchanged. The trace reports it
#' as \code{"hessian modified (capped)"}. The same scaling is applied to the
#' gradient the method falls back on when a solve fails.
#'
#' Without it the length of a repaired step is unbounded and the line search is
#' the only thing that bounds it, at one objective evaluation per backtrack:
#' measured on a marginal criterion whose outer Hessian is indefinite at
#' ordinary points, a gradient of 29 along a direction floored at
#' \eqn{10^{-8}\lambda_{\max}} gave a step of \eqn{4\times 10^{4}} on a log
#' scale, and each of the 23 backtracks that followed was a whole penalized
#' refit that could not be evaluated.
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
#' @references
#' Gill, P. E., Murray, W. and Wright, M. H. (1981).
#' \emph{Practical Optimization}. Academic Press, London.
#'
#' Nocedal, J. and Wright, S. J. (2006). \emph{Numerical Optimization},
#' 2nd edition. Springer, New York.
#'
#' @export
newton <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
                   hessian_mod = c("eigen", "ridge"), floor = 1e-8,
                   step = 1, line_search = armijo(),
                   maxit = 200, max_eval = Inf,
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
#' @references
#' Broyden, C. G. (1970). The convergence of a class of double-rank
#' minimization algorithms. \emph{IMA Journal of Applied Mathematics}
#' \strong{6}, 76--90. The update was obtained independently the same
#' year by Fletcher, Goldfarb and Shanno, whence the name.
#'
#' Nocedal, J. and Wright, S. J. (2006). \emph{Numerical Optimization},
#' 2nd edition. Springer, New York.
#'
#' @export
bfgs <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
                 curv_tol = 1e-10, max_skip = 5,
                 step = 1, line_search = wolfe(),
                 maxit = 500, max_eval = Inf,
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
#' @references
#' Nocedal, J. (1980). Updating quasi-Newton matrices with limited
#' storage. \emph{Mathematics of Computation} \strong{35}, 773--782.
#'
#' Liu, D. C. and Nocedal, J. (1989). On the limited memory BFGS method
#' for large scale optimization. \emph{Mathematical Programming}
#' \strong{45}, 503--528.
#'
#' @export
lbfgs <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
                  memory = 10, curv_tol = 1e-10,
                  step = 1, line_search = wolfe(),
                  maxit = 500, max_eval = Inf,
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

#' @title Minimize by Newton's Method
#' @name minimize.Newton
#' @param optimizer A \code{Newton} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Newton) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    # A numerical Hessian differences the gradient once per coordinate and
    # direction; without an analytic gradient each of those is itself 2p
    # objective values, so one iteration costs about 4p^2 evaluations. When the
    # evaluation budget admits fewer than two such iterations the run can only
    # end on the budget, so that is said here rather than discovered from a
    # one-iteration result.
    if (is.null(he) && is.function(fn)) {
      p <- length(par)
      per_iter <- if (is.null(gr)) 4 * p^2 + 2 * p else 0
      if (per_iter > 0 && 2 * per_iter > optimizer@max_eval) {
        warning("newton() without 'gr' and 'he' costs about 4*p^2 = ",
                format(4 * p^2), " objective evaluations per iteration at p = ",
                p, ",
  and max_eval = ", format(optimizer@max_eval),
                " admits fewer than two iterations. Supply 'gr' (or 'he'),
",
                "  raise 'max_eval', or use bfgs(), which needs no Hessian.",
                call. = FALSE)
      }
    }
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "newton", hessian_mod = optimizer@hessian_mod,
                     floor = optimizer@floor))
  }

#' @title Minimize by BFGS
#' @name minimize.Bfgs
#' @param optimizer A \code{Bfgs} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Bfgs) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "bfgs", curv_tol = optimizer@curv_tol,
                     max_skip = as.integer(optimizer@max_skip)))
  }

#' @title Minimize by Limited-Memory BFGS
#' @name minimize.Lbfgs
#' @param optimizer An \code{Lbfgs} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Lbfgs) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "lbfgs", memory = as.integer(optimizer@memory),
                     curv_tol = optimizer@curv_tol))
  }
