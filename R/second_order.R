#' @include optimizer_class.R
#' @include generics.R
#' @include line_search.R
#' @include modify.R
NULL

# Newton, BFGS and L-BFGS. Each is a constructor and a method that names its
# direction; the loop, the line search, the stopping rule and the reporting are
# all shared, which is what the frame was built for.


#' @title S7 Class for Newton's Method
#'
#' @description
#' This class is an optimizer that holds the initial step length, the line
#' search, the [hessian_modification()] that repairs an indefinite Hessian,
#' and the two settings that bound and scale a step. [newton()] builds it. Among the shipped methods, only this one and [chebyshev()], which inherits from it,
#' read the `he` argument of [minimize()].
#'
#' @details
#' Beyond the seven properties of every optimizer, a `Newton` carries five:
#' `step` and `line_search`, shared with the other line-search methods, and
#' `modify`, `max_length` and `typical`, which are its own.
#'
#' @param modify The [hessian_modification()] object that repairs an
#'   indefinite Hessian; see [newton()].
#' @param max_length The largest component that a step may have.
#' @param typical `NULL`, or the typical size of each parameter, in the
#'   units of which the step lengths are measured.
#' @param step,line_search The initial step length and the [line_search()]
#'   object, as in [newton()].
#'
#' @return An S7 object of class `Newton` inheriting from [optimizer()], with
#'   the five properties above beside the seven shared ones.
#'
#' @seealso [newton()] for the constructor, [bfgs()] for the method that
#'   needs no Hessian.
#' @name Newton-class
#' @aliases Newton
#' @keywords internal
Newton <- S7::new_class("Newton", parent = optimizer,
  properties = list(step = S7::class_numeric, line_search = S7::class_any,
                    modify = S7::class_any,
                    max_length = S7::class_numeric,
                    typical = S7::class_any))


#' @title Newton's Method with a Modified Hessian
#'
#' @description
#' Solves \eqn{H d = -g} for the direction, repairing \eqn{H} when it is not
#' positive definite, and then performs a line search along it. The method
#' converges quadratically near a minimum.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`; see
#'   [crit_any()] for the effect of that disjunction.
#' @param modify A [hessian_modification()] object that determines how an
#'   indefinite or nearly singular Hessian is repaired before the step. A
#'   character string is rejected. The constructors are:
#'   * [modify_eigen()] (default): raises the eigenvalues of the Hessian that
#'     are below a floor to the floor.
#'   * [modify_ridge()]: adds a multiple of the identity, doubled until the
#'     Cholesky factorization succeeds.
#'   * A user-defined subclass of [hessian_modification()] with a method for
#'     [hessian_modify()], which needs no change to this function.
#'
#'   The floor is an argument of the constructor, as in
#'   `newton(modify = modify_ridge(floor = 1e-6))`. See [modifications] for
#'   the comparison and [hessian_modify()] for writing a new object.
#' @param step Initial step length offered to the line search, a single
#'   positive number. Defaults to `1`, which is the Newton step itself and is
#'   accepted unchanged near the solution.
#' @param line_search A [line_search()] object. Defaults to [armijo()]. The
#'   curvature condition that [bfgs()] needs is unnecessary here, because the
#'   curvature is read from the Hessian.
#' @param maxit,max_eval,verbose,refresh,keep_trace As in [optimizer()].
#'   `maxit` defaults to 200 here and to 500 in [bfgs()] and [lbfgs()]: a
#'   Newton run that has not converged in 200 iterations usually has a
#'   problem that more iterations do not solve.
#' @param max_length The largest component that a step may have, a single
#'   positive number. A direction longer than this in the infinity norm is
#'   scaled down to this length before the line search. Defaults to `Inf`,
#'   which leaves every step as it is. See below.
#' @param typical `NULL` (the default), or the typical size of each
#'   parameter: one positive number, or one per parameter. Every length used
#'   by the method is then measured in these units, as the largest
#'   \eqn{\lvert d_i\rvert / t_i} in place of the largest
#'   \eqn{\lvert d_i\rvert}: the bound `max_length`, and the scaling of a
#'   direction whose length is not set by the curvature. With `NULL` the
#'   lengths are measured in the units of the parameters.
#'
#' @details
#' # Why the Hessian has to be repaired
#'
#' \eqn{-H^{-1}g} is guaranteed to be a descent direction only when \eqn{H}
#' is positive definite. Where the objective curves downwards the unmodified
#' step can point towards a saddle or a maximum, and a line search along an
#' ascent direction cannot decrease the objective. Indefinite Hessians are
#' common far from the solution.
#'
#' The repair is selected by `modify`. Every repair begins with an attempted
#' Cholesky factorization, which, when it succeeds, is both the test for
#' positive definiteness and the solve. When it fails, the Hessian is passed
#' to [hessian_modify()] and the two shipped repairs are:
#'
#' \describe{
#'   \item{[modify_eigen()]}{decompose \eqn{H} and raise every eigenvalue below
#'     \eqn{\max(\texttt{floor}, \texttt{floor}\,\max_i\lvert\lambda_i\rvert)}
#'     to that value. The direction is then the Newton direction in the
#'     subspace where the curvature is positive, and close to the gradient
#'     direction in the rest. The repair costs a symmetric eigendecomposition
#'     and gives the best-conditioned matrix. If the eigendecomposition
#'     fails, the ridge repair is used.}
#'   \item{[modify_ridge()]}{add \eqn{\tau I} and double \eqn{\tau} until the
#'     factorization succeeds. This is Levenberg's idea: it interpolates
#'     between the Newton step at \eqn{\tau = 0} and a scaled steepest-descent
#'     step for large \eqn{\tau}. It costs less than the eigendecomposition and
#'     shifts every eigenvalue, not only the small ones.}
#' }
#'
#' A user-defined [hessian_modification()] is called through R at each such
#' iterate. A repair of the Hessian is recorded in the trace as
#' `hessian modified`, or `hessian modified (capped)` when the repaired step
#' is also scaled down (see below). The label does not record which repair
#' was used.
#'
#' # Bounding the length of a step
#'
#' The Newton step minimizes the local quadratic model of the objective, and
#' its length is set by the curvature at the current point. Where the
#' objective flattens towards an asymptote the curvature is small and the
#' step is long, and it can carry the iterate past the minimum into the flat
#' region. The line search accepts such a step whenever it lowers the
#' objective, and in the flat region the gradient is close to zero, so the
#' run can stop there. Such steps occur, for example, when a smoothing
#' parameter is optimized on the logarithmic scale.
#'
#' `max_length` bounds the step; it is the `maxstep` of the line-search
#' Newton method of Dennis and Schnabel (1983) and the `maxNstep` of the
#' smoothing-parameter iteration of mgcv. The bound only shortens a step, and
#' near a minimum the Newton steps are short, so the quadratic convergence of
#' the method is kept. A step shortened by the bound is recorded in the trace
#' as `newton step capped`, except at an iteration at which the Hessian was
#' repaired, which keeps the label of the repair.
#'
#' # Parameters in different units
#'
#' A length in the units of the parameters has different meanings for
#' parameters on different scales: a displacement of one is a long step for
#' a parameter on the log scale and a short one for a variance of several
#' hundred. `typical` gives the typical size \eqn{t_i} of each parameter, as
#' in the scaled norm of Dennis and Schnabel (1983): every length is then the
#' largest \eqn{\lvert d_i\rvert / t_i}, for the bound `max_length` and for
#' the scaling of a repaired direction or of a gradient fallback.
#'
#' # The length of a repaired step
#'
#' Where the factorization succeeds, \eqn{-H^{-1}g} is the Newton step and
#' `step = 1` is its natural unit, so the step is taken as it stands, subject
#' only to `max_length`. Where it fails, the component of the direction in
#' the floored subspace is \eqn{g_i/\lambda_{\mathrm{floor}}}, whose size is
#' set by the `floor` of the modification and not by the curvature of the
#' objective. The direction is
#' therefore multiplied by \eqn{\min(1, 1/\lVert d\rVert_\infty)}, which
#' gives a displacement of order one in the parameters. [bfgs()] and
#' [lbfgs()] apply the same rule to their first direction. The factor never
#' exceeds one, so a run whose repaired steps were already of order one is
#' unchanged. The same scaling applies to the gradient direction that the
#' method uses when a solve fails.
#'
#' Without this factor the length of a repaired step is bounded only by the
#' line search, at one objective evaluation per backtrack. Each backtrack can
#' be expensive, for example a complete penalized refit when the objective is
#' a marginal criterion.
#'
#' # The Hessian by differences
#'
#' When `he` is not supplied to [minimize()] the Hessian is obtained by
#' differencing the gradient. With an analytic gradient that is one numerical
#' differentiation. When the gradient is *also* differenced, the Hessian is a
#' difference of differences, which happens only in [newton()] and
#' [chebyshev()]: one Newton iteration then costs about \eqn{4p^2}
#' evaluations of the objective, against about \eqn{2p} evaluations of the
#' gradient with an analytic gradient. For an objective without either
#' derivative, [bfgs()], which needs no Hessian, is usually cheaper.
#'
#' @return An S7 object of class [Newton], inheriting from [optimizer()], to
#'   be handed to [minimize()].
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
#' # The repair of an indefinite Hessian is an object, with its floor.
#' newton(modify = modify_ridge(floor = 1e-6))
#'
#' # At a saddle the unrepaired step points the wrong way. The Hessian at
#' # the start has eigenvalues 2 and -1.99; the repair is recorded in the
#' # trace.
#' sad <- function(p) p[1]^2 - p[2]^2 + 0.1 * p[2]^4
#' sg  <- function(p) c(2 * p[1], -2 * p[2] + 0.4 * p[2]^3)
#' sh  <- function(p) matrix(c(2, 0, 0, -2 + 1.2 * p[2]^2), 2, 2)
#' eigen(sh(c(0.5, 0.1)), only.values = TRUE)$values
#'
#' r <- minimize(newton(keep_trace = TRUE), sad, c(0.5, 0.1), gr = sg, he = sh)
#' c(value = r@value, iterations = r@iterations)
#' unique(r@trace$safeguard)
#'
#' # summary() counts the safeguards that fired during the run.
#' summary(minimize(newton(keep_trace = TRUE), rosen, c(-1.2, 1),
#'                  gr = rosen_gr, he = rosen_he))
#'
#' @seealso [modifications] for the repairs of the Hessian, [modify_eigen()],
#'   [modify_ridge()] and [hessian_modify()] to write another,
#'   [bfgs()] and [lbfgs()] when a Hessian is not available,
#'   [armijo()] for the line search, [summary.optimizer_result()] for the
#'   safeguard counts.
#' @references
#' Dennis, J. E. and Schnabel, R. B. (1983). *Numerical Methods for
#' Unconstrained Optimization and Nonlinear Equations*. Prentice-Hall,
#' Englewood Cliffs, NJ.
#'
#' Gill, P. E., Murray, W. and Wright, M. H. (1981).
#' *Practical Optimization*. Academic Press, London.
#'
#' Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*,
#' 2nd edition. Springer, New York.
#'
#' @export
newton <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
                   modify = modify_eigen(),
                   step = 1, line_search = armijo(),
                   maxit = 200, max_eval = Inf,
                   verbose = FALSE, refresh = 10, keep_trace = FALSE,
                   max_length = Inf, typical = NULL) {
  check_modify(modify)
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  if (length(max_length) != 1L || !is.numeric(max_length) ||
      is.na(max_length) || max_length <= 0) {
    stop("'max_length' must be a single positive number (Inf for no bound).",
         call. = FALSE)
  }
  if (!is.null(typical) &&
      (!length(typical) || !is.numeric(typical) || anyNA(typical) ||
       any(!is.finite(typical)) || any(typical <= 0))) {
    stop(paste0("'typical' must be NULL, or positive finite numbers: one, ",
                "or one per parameter."), call. = FALSE)
  }
  Newton(name = "Newton", criterion = criterion, maxit = maxit,
         max_eval = max_eval, verbose = verbose, refresh = refresh,
         keep_trace = keep_trace, step = step, line_search = line_search,
         modify = modify, max_length = max_length,
         typical = typical)
}


#' @title S7 Class for BFGS
#'
#' @description
#' This class is an optimizer that holds the initial step length, the line
#' search, and the two settings that apply when a secant pair has no usable
#' curvature. [bfgs()] builds it. The inverse-Hessian approximation is not a
#' property: it is built inside the run and discarded at its end.
#'
#' @details
#' Beyond the seven properties of every optimizer, a `Bfgs` carries four:
#' `step` and `line_search`, shared with the other line-search methods, and
#' `curv_tol` and `max_skip`, which are its own.
#'
#' @param curv_tol The curvature threshold below which the update is skipped.
#' @param max_skip Consecutive skips before the approximation is reset to the
#'   identity.
#' @param step,line_search The initial step length and the [line_search()]
#'   object, as in [bfgs()].
#'
#' @return An S7 object of class `Bfgs` inheriting from [optimizer()], with
#'   the four properties above beside the seven shared ones.
#'
#' @seealso [bfgs()] for the constructor, [Lbfgs] for the limited-memory
#'   version.
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
#' Builds an approximation to the inverse Hessian from successive gradients,
#' so the direction is a matrix-vector product and the method never
#' evaluates second derivatives. It is the usual choice for a smooth problem of moderate
#' size.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`.
#' @param curv_tol The update is skipped when
#'   \eqn{s^\top y \le \texttt{curv\_tol}\,\lVert s\rVert\lVert y\rVert}.
#'   Defaults to `1e-10`. A larger value skips more often. Zero is admitted
#'   and gives the textbook condition \eqn{s^\top y > 0}. A negative value is
#'   rejected, because it would admit pairs with negative curvature, and the
#'   approximation would then no longer be guaranteed to stay positive
#'   definite.
#' @param max_skip Consecutive skipped updates before the approximation is
#'   reset to the identity. Defaults to 5, and has to be a whole number.
#'   Zero is admitted and resets on the first skip; a negative value is
#'   rejected, because the comparison `>=` would treat it as zero.
#' @param step,line_search,maxit,max_eval,verbose,refresh,keep_trace As in
#'   [newton()], except that `line_search` defaults to [wolfe()] and `maxit`
#'   to 500. See below for the reason for the line search.
#'
#' @details
#' # Notation
#'
#' \eqn{s = x_{new} - x_{old}} is the secant vector and
#' \eqn{y = g_{new} - g_{old}} the change in the gradient. The scalar step
#' length that a line search returns is \eqn{\alpha}, a different quantity;
#' the two are related by \eqn{s = \alpha d} for the direction \eqn{d}.
#'
#' # Why the line search is the strong Wolfe one
#'
#' The update is meaningful only when \eqn{s^\top y > 0}, and the Wolfe
#' curvature condition guarantees that it holds: Armijo backtracking can
#' accept a step so short that the gradient has barely changed, which gives
#' a pair without curvature.
#'
#' On the problems of [test_problems()] the difference is small in practice.
#' Under [armijo()], BFGS skips an update on a few problems, needs a few more
#' iterations on two of them, and converges on the non-smooth `abs_sum`,
#' where it stops without converging under [wolfe()]. Wolfe is the default
#' because it guarantees a sound update; Armijo is a reasonable choice where
#' evaluations are expensive.
#'
#' # Skipping and resetting
#'
#' When the curvature condition fails, the update is skipped. A small
#' \eqn{s^\top y} makes \eqn{\rho = 1/s^\top y} very large, and one such
#' update can destroy the accumulated approximation; skipping it keeps the
#' previous matrix, which is still positive definite. After `max_skip`
#' consecutive skips the matrix is reset to the identity, so the method
#' restarts as steepest descent and rebuilds the approximation. Both events
#' appear in the trace, as `bfgs update skipped` and `bfgs reset`. After a
#' reset of this kind the next direction is \eqn{-g}, without the scaling
#' described below; after a reset caused by a direction that is not a
#' descent direction, the direction is scaled.
#'
#' The first accepted pair rescales the identity by \eqn{s^\top y / y^\top y}.
#' Without it the first quasi-Newton step is taken with a unit Hessian, which
#' on a badly scaled problem has the wrong magnitude and costs backtracks.
#' Until a pair has been accepted, the direction is the negative gradient
#' multiplied by \eqn{\min(1, 1/\lVert g\rVert_\infty)}, which gives a
#' displacement of order one in the parameters and never lengthens the
#' step.
#'
#' @return An S7 object of class [Bfgs], inheriting from [optimizer()], to be
#'   handed to [minimize()].
#'
#' @examples
#' rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' rg <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                     200 * (p[2] - p[1]^2))
#' minimize(bfgs(), rosen, c(-1.2, 1), gr = rg)
#'
#' # With Armijo instead of Wolfe, the update is skipped once here, the trace
#' # records it, and the run converges.
#' a <- minimize(bfgs(line_search = armijo(), keep_trace = TRUE), rosen,
#'               c(-1.2, 1), gr = rg)
#' c(iterations = a@iterations, skips = sum(a@trace$safeguard == "bfgs update skipped"))
#'
#' # A curvature threshold that every pair fails to meet forces every update
#' # to be skipped; the run exhausts its budget, and the trace names both
#' # safeguards.
#' b <- minimize(bfgs(curv_tol = 1e10, keep_trace = TRUE), rosen, c(-1.2, 1),
#'               gr = rg)
#' c(converged = b@converged, value = b@value)
#' unique(b@trace$safeguard)
#'
#' @seealso [lbfgs()] for many parameters, [newton()] when a Hessian is
#'   available, [wolfe()] for the default line search of this method.
#' @references
#' Broyden, C. G. (1970). The convergence of a class of double-rank
#' minimization algorithms. *IMA Journal of Applied Mathematics*
#' **6**, 76--90. The update was obtained independently the same
#' year by Fletcher, Goldfarb and Shanno, whence the name.
#'
#' Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*,
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
  check_nonneg(curv_tol, "curv_tol")
  check_whole(max_skip, "max_skip")
  Bfgs(name = "BFGS", criterion = criterion, maxit = maxit,
       max_eval = max_eval, verbose = verbose, refresh = refresh,
       keep_trace = keep_trace, step = step, line_search = line_search,
       curv_tol = curv_tol, max_skip = max_skip)
}


#' @title S7 Class for Limited-Memory BFGS
#'
#' @description
#' This class is an optimizer that holds the number of secant pairs to keep,
#' the curvature threshold that a pair must meet to be stored, and the step
#' length and the line search. [lbfgs()] builds it. The pairs themselves
#' exist only inside the run.
#'
#' @details
#' Beyond the seven properties of every optimizer, an `Lbfgs` carries four:
#' `step` and `line_search`, and `memory` and `curv_tol`, which are its own.
#' It has no counterpart of `max_skip` of [Bfgs]: a pair that fails the
#' curvature test is not stored (the trace records `lbfgs update skipped`),
#' and there is no accumulated matrix to protect. When the direction is not a
#' descent direction, the stored pairs are cleared (`lbfgs reset`).
#'
#' @param memory How many secant pairs to keep.
#' @param curv_tol The curvature threshold below which a pair is not stored.
#' @param step,line_search The initial step length and the [line_search()]
#'   object, as in [lbfgs()].
#'
#' @return An S7 object of class `Lbfgs` inheriting from [optimizer()], with
#'   the four properties above beside the seven shared ones.
#'
#' @seealso [lbfgs()] for the constructor, [Bfgs] for the full-matrix
#'   version.
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
#' Applies the BFGS update without forming the matrix: only the last `memory`
#' secant pairs are kept, and the direction comes from the two-loop
#' recursion. It costs \eqn{O(mp)} in time and memory where [bfgs()] costs
#' \eqn{O(p^2)}. The
#' direction equals the one of the explicitly assembled inverse built from
#' the same pairs and the same initial scaling. The two methods do not give
#' the same directions in a run, because [bfgs()] keeps every pair and scales
#' its initial matrix by the first accepted pair, whereas the recursion uses
#' the last `memory` pairs and the most recent one for the scaling.
#'
#' @param criterion The stopping rule, a [criterion()] object. Defaults to
#'   `crit_any(crit_grad(), crit_abs_obj(), crit_abs_par())`.
#' @param memory How many secant pairs to keep, a single positive whole
#'   number. Defaults to 10. A fractional value is rejected instead of
#'   truncated, since the count is passed to the compiled code as an integer.
#' @param curv_tol A pair is not stored when
#'   \eqn{s^\top y \le \texttt{curv\_tol}\,\lVert s\rVert\lVert y\rVert}.
#'   Defaults to `1e-10`. Zero is admitted and gives \eqn{s^\top y > 0}; a
#'   negative value is rejected, for the reason given on the page of
#'   [bfgs()].
#' @param step,line_search,maxit,max_eval,verbose,refresh,keep_trace As in
#'   [bfgs()].
#'
#' @details
#' # Notation
#'
#' \eqn{s = x_{new} - x_{old}} is the secant vector and
#' \eqn{y = g_{new} - g_{old}} the change in the gradient. The scalar step
#' length is \eqn{\alpha}, a different quantity.
#'
#' # The saving
#'
#' The two-loop recursion returns exactly the product of the explicitly
#' assembled inverse built from the same pairs and the same scaling (the
#' last example). The methods differ in cost: one iteration costs
#' \eqn{O(p^2)} for [bfgs()] and \eqn{O(mp)} for [lbfgs()], while the
#' numbers of iterations of the two methods on a dense quadratic are close.
#' The difference in time becomes large for a few hundred parameters, and
#' for several thousand the \eqn{O(p^2)} storage and work of the full matrix
#' become prohibitive.
#'
#' The recursion is scaled at each iteration by
#' \eqn{s^\top y / y^\top y} of the most recent pair. This scalar estimates
#' the inverse curvature along the latest step and takes the place of the
#' initial matrix of the full update.
#'
#' # Choosing memory
#'
#' A larger `memory` is not always better. Old pairs describe the curvature
#' at points that the iterate has left, and each pair costs \eqn{O(p)} per
#' iteration. Beyond a handful of pairs the number of iterations changes
#' little while the arithmetic keeps rising; ten is the conventional
#' choice.
#'
#' @return An S7 object of class [Lbfgs], inheriting from [optimizer()], to
#'   be handed to [minimize()].
#'
#' @examples
#' rosen <- function(p) (1 - p[1])^2 + 100 * (p[2] - p[1]^2)^2
#' rg <- function(p) c(-2 * (1 - p[1]) - 400 * p[1] * (p[2] - p[1]^2),
#'                     200 * (p[2] - p[1]^2))
#' minimize(lbfgs(memory = 5), rosen, c(-1.2, 1), gr = rg)
#'
#' # The two-loop recursion is not an approximation to the full update: on
#' # the same pairs it returns the same product, to rounding error.
#' set.seed(1)
#' p <- 6
#' S <- matrix(rnorm(p * 4), p, 4)
#' Y <- matrix(rnorm(p * 4), p, 4) + 3 * S      # positive curvature
#' g <- rnorm(p)
#' gamma <- sum(S[, 4] * Y[, 4]) / sum(Y[, 4]^2)
#'
#' H <- diag(gamma, p)                          # the explicit inverse
#' for (j in 1:4) {
#'   s <- S[, j]; y <- Y[, j]; rho <- 1 / sum(s * y)
#'   V <- diag(p) - rho * outer(s, y)
#'   H <- V %*% H %*% t(V) + rho * outer(s, s)
#' }
#'
#' two_loop <- function(g, S, Y, gamma) {       # the recursion
#'   m <- ncol(S); a <- numeric(m); q <- g
#'   for (j in m:1) {
#'     rho <- 1 / sum(S[, j] * Y[, j])
#'     a[j] <- rho * sum(S[, j] * q); q <- q - a[j] * Y[, j]
#'   }
#'   r <- gamma * q
#'   for (j in 1:m) {
#'     rho <- 1 / sum(S[, j] * Y[, j])
#'     b <- rho * sum(Y[, j] * r); r <- r + S[, j] * (a[j] - b)
#'   }
#'   r
#' }
#' max(abs(H %*% g - two_loop(g, S, Y, gamma)))
#'
#' @seealso [bfgs()] for a few parameters, [newton()] when a Hessian is
#'   available.
#' @references
#' Nocedal, J. (1980). Updating quasi-Newton matrices with limited
#' storage. *Mathematics of Computation* **35**, 773--782.
#'
#' Liu, D. C. and Nocedal, J. (1989). On the limited memory BFGS method
#' for large scale optimization. *Mathematical Programming*
#' **45**, 503--528.
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
  check_nonneg(curv_tol, "curv_tol")
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
#' @param optimizer A `Newton` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()].
#' @return An [optimizer_result()].
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
        warning("newton() without 'gr' and 'he' costs about 4*p^2 + 2*p = ",
                format(per_iter), " objective evaluations per iteration at p = ",
                p, ",
  and max_eval = ", format(optimizer@max_eval),
                " admits fewer than two iterations. Supply 'gr' (or 'he'),
",
                "  raise 'max_eval', or use bfgs(), which needs no Hessian.",
                call. = FALSE)
      }
    }
    typ <- optimizer@typical
    if (!is.null(typ)) {
      if (length(typ) == 1L) typ <- rep(typ, length(par))
      if (length(typ) != length(par)) {
        stop(sprintf(paste0("'typical' has %d sizes and 'par' has %d ",
                            "values; give one size, or one per parameter."),
                     length(typ), length(par)), call. = FALSE)
      }
    }
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                c(list(type = "newton"), modify_spec(optimizer@modify),
                  list(max_length = optimizer@max_length,
                       typical = if (is.null(typ)) numeric(0) else typ)))
  }

#' @title Minimize by BFGS
#' @name minimize.Bfgs
#' @param optimizer A `Bfgs` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()].
#' @return An [optimizer_result()].
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
#' @param optimizer An `Lbfgs` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()].
#' @return An [optimizer_result()].
#' @keywords internal
S7::method(minimize, Lbfgs) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "lbfgs", memory = as.integer(optimizer@memory),
                     curv_tol = optimizer@curv_tol))
  }
