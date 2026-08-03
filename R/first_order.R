#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

# The methods that use the gradient and nothing else. All three are a direction
# plus a line search, so they share the compiled loop of descent.cpp and differ
# only in how the direction is formed.


#' @title S7 Class for Gradient Descent
#' @description The class \code{\link{gd}} instantiates.
#' @param step The initial step length offered to the line search.
#' @param line_search A \code{\link{line_search}} object.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{gd}}
#' @name GradientDescent-class
#' @aliases GradientDescent
#' @keywords internal
GradientDescent <- S7::new_class("GradientDescent", parent = optimizer,
  properties = list(
    step        = S7::class_numeric,
    line_search = S7::class_any
  ))


#' @title Gradient Descent
#'
#' @description
#' The simplest method there is: step along the negative gradient, and let the
#' line search decide how far.
#'
#' @param criterion The stopping rule; see \code{\link{crit_any}}.
#' @param step Initial step length offered to the line search each iteration.
#'   Defaults to \code{1}.
#' @param line_search How far to go along the direction; see
#'   \code{\link{armijo}} and \code{\link{wolfe}}. Defaults to Armijo
#'   backtracking, which is all a method with nothing to protect needs.
#' @param maxit Maximum iterations. Defaults to 500.
#' @param max_eval Maximum objective evaluations. Defaults to 10000.
#' @param verbose Report progress? Defaults to \code{FALSE}.
#' @param refresh Report every this many iterations. Defaults to 10.
#' @param keep_trace Store the iteration path? Defaults to \code{FALSE}.
#'
#' @details
#' The direction \eqn{-g} minimises \eqn{g^\top d} over directions of a given
#' Euclidean length. Under an exact line search consecutive directions are
#' orthogonal, so on an ill-conditioned objective the iterates zigzag and
#' converge slowly; \code{\link{cg}} corrects this by combining each direction
#' with the previous one at no extra cost per iteration. Gradient descent is
#' mainly useful as a baseline: on smooth problems \code{\link{bfgs}} converges
#' in far fewer iterations.
#'
#' @return An S7 object of class \code{GradientDescent}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @examples
#' gd()
#' gd(criterion = crit_grad(1e-10), maxit = 2000)
#'
#' minimize(gd(), function(p) sum((p - c(1, 2))^2), c(0, 0),
#'          gr = function(p) 2 * (p - c(1, 2)))
#'
#' @seealso \code{\link{cg}}, \code{\link{bb}}, \code{\link{bfgs}}
#' @export
gd <- function(criterion = crit_any(crit_grad(1e-8), crit_rel_obj(1e-12)),
               step = 1, line_search = armijo(),
               maxit = 500, max_eval = 10000,
               verbose = FALSE, refresh = 10, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  GradientDescent(
    name = "gradient descent", criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    step = step, line_search = line_search
  )
}


#' @title Minimise by Gradient Descent
#' @name minimize.GradientDescent
#' @description Runs \code{\link{gd}} on the objective.
#' @param optimizer A \code{GradientDescent} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, GradientDescent) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    run_descent(optimizer, fn, par, gr, he, lower, upper, list(type = "gd"))
  }


# --- conjugate gradients ----------------------------------------------------

#' @title S7 Class for Conjugate Gradients
#' @description The class \code{\link{cg}} instantiates.
#' @param beta Which update formula.
#' @param restart_every How often the method is restarted at steepest descent.
#' @param step The initial step length offered to the line search.
#' @param line_search A \code{\link{line_search}} object.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{cg}}
#' @name Cg-class
#' @aliases Cg
#' @keywords internal
Cg <- S7::new_class("Cg", parent = optimizer,
  properties = list(
    beta          = S7::class_character,
    restart_every = S7::class_numeric,
    step          = S7::class_numeric,
    line_search   = S7::class_any
  ))


#' @title Nonlinear Conjugate Gradients
#'
#' @description
#' Gradient descent with the zigzag taken out: each direction is bent back
#' towards the previous one. It stores two vectors and no matrix, which makes it
#' the classical answer for a problem too large to hold a Hessian.
#'
#' @param criterion The stopping rule; see \code{\link{crit_any}}.
#' @param beta Which formula for the bend: \code{"pr"} (Polak–Ribière, the
#'   default), \code{"fr"} (Fletcher–Reeves), \code{"hs"} (Hestenes–Stiefel) or
#'   \code{"dy"} (Dai–Yuan). See Details.
#' @param restart_every Restart at steepest descent every this many iterations.
#'   Defaults to \code{0}, meaning never; a positive value is usually the
#'   dimension of the problem.
#' @param step Initial step length offered to the line search. Defaults to
#'   \code{1}.
#' @param line_search Defaults to \code{wolfe(c2 = 0.1)}; the convergence
#'   theory of the method assumes a strong Wolfe step, see Details.
#' @param maxit Maximum iterations. Defaults to 1000.
#' @param max_eval Maximum objective evaluations. Defaults to 20000.
#' @param verbose Report progress? Defaults to \code{FALSE}.
#' @param refresh Report every this many iterations. Defaults to 20.
#' @param keep_trace Store the iteration path? Defaults to \code{FALSE}.
#'
#' @details
#' The direction is
#' \deqn{d_k = -g_k + \beta_k d_{k-1},}
#' and the whole method is the choice of \eqn{\beta}. On a quadratic with an
#' exact line search the directions come out conjugate with respect to the
#' Hessian, so the method terminates in \eqn{p} steps exactly — \emph{without
#' ever forming that Hessian}, which is the point. The storage is two vectors
#' against \code{\link{bfgs}}'s \eqn{p \times p} matrix.
#'
#' \subsection{Choice of beta}{
#' All four agree on a quadratic with an exact line search and differ everywhere
#' else. \code{"fr"} has the cleanest convergence theory and the well-known
#' practical fault of stalling for many iterations after a poor step.
#' \code{"pr"} recovers from a poor step immediately, because a small
#' \eqn{y = g_k - g_{k-1}} sends \eqn{\beta} towards zero and the method back to
#' steepest descent; its known theoretical non-convergence is repaired by
#' clamping \eqn{\beta} at zero, which restarts the method and is recorded as a
#' restart in the trace. \code{"hs"} and \code{"dy"} are the other two standard
#' choices.
#' }
#'
#' \subsection{Line search}{
#' The theory behind every one of these formulas assumes a step satisfying the
#' \strong{strong} Wolfe conditions, and uses it to prove that the direction
#' produced is a descent direction at all. Backtracking gives no such guarantee,
#' so \code{\link{wolfe}} is the default and departing from it is departing
#' from the theory.
#'
#' The constant matters as much as the search. \code{\link{wolfe}} defaults to
#' \eqn{c_2 = 0.9}, which is right for \code{\link{bfgs}}, where the curvature
#' approximation repairs a loose step. Conjugate gradients has nothing to repair
#' with: the accumulated conjugacy is only as good as the line search that
#' produced it, and a loose one degrades it directly. Measured on Rosenbrock,
#' \eqn{c_2 = 0.9} needs 120 iterations against 35 at \eqn{c_2 = 0.1}, which is
#' why the default here is the tighter one.
#'
#' As a safeguard against the cases the theory misses, a direction that comes
#' out non-descent is replaced by \eqn{-g} and the substitution is reported.
#' }
#'
#' @return An S7 object of class \code{Cg}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @references
#' Hestenes, M. R. and Stiefel, E. (1952). Methods of conjugate gradients for
#' solving linear systems. \emph{Journal of Research of the NBS} \strong{49},
#' 409--436.
#'
#' Polak, E. and Ribière, G. (1969). Note sur la convergence de méthodes de
#' directions conjuguées. \emph{Revue Française d'Informatique et de Recherche
#' Opérationnelle} \strong{3}, 35--43.
#'
#' Dai, Y. H. and Yuan, Y. (1999). A nonlinear conjugate gradient method with a
#' strong global convergence property. \emph{SIAM Journal on Optimization}
#' \strong{10}, 177--182.
#'
#' @examples
#' cg()
#' cg(beta = "fr", restart_every = 10)
#'
#' f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
#' gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
#'                     200 * (p[2] - p[1]^2))
#' minimize(cg(), f, c(-1.2, 1), gr = gr)@par
#'
#' @seealso \code{\link{gd}}, \code{\link{lbfgs}}, \code{\link{bb}}
#' @export
cg <- function(criterion = crit_any(crit_grad(1e-8), crit_rel_obj(1e-12)),
               beta = c("pr", "fr", "hs", "dy"), restart_every = 0,
               step = 1, line_search = wolfe(c2 = 0.1),
               maxit = 1000, max_eval = 20000,
               verbose = FALSE, refresh = 20, keep_trace = FALSE) {
  beta <- match.arg(beta)
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  if (length(restart_every) != 1L || !is.numeric(restart_every) ||
      is.na(restart_every) || restart_every < 0) {
    stop("'restart_every' must be a single non-negative number.", call. = FALSE)
  }
  Cg(
    name = paste0("conjugate gradients (", beta, ")"), criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    beta = beta, restart_every = restart_every,
    step = step, line_search = line_search
  )
}


#' @title Minimise by Conjugate Gradients
#' @name minimize.Cg
#' @description Runs \code{\link{cg}} on the objective.
#' @param optimizer A \code{Cg} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Cg) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "cg", beta = optimizer@beta,
                     restart_every = as.integer(optimizer@restart_every)))
  }


# --- Barzilai-Borwein -------------------------------------------------------

#' @title S7 Class for the Barzilai-Borwein Method
#' @description The class \code{\link{bb}} instantiates.
#' @param variant Which step-length formula.
#' @param alpha0 The step length used before there is a secant pair.
#' @param alpha_min,alpha_max Bounds on it.
#' @param curv_tol The relative threshold below which a pair is refused.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{bb}}
#' @name Bb-class
#' @aliases Bb
#' @keywords internal
Bb <- S7::new_class("Bb", parent = optimizer,
  properties = list(
    variant     = S7::class_character,
    alpha0      = S7::class_numeric,
    alpha_min   = S7::class_numeric,
    alpha_max   = S7::class_numeric,
    curv_tol    = S7::class_numeric,
    step        = S7::class_numeric,
    line_search = S7::class_any
  ))


#' @title The Barzilai-Borwein Method
#'
#' @description
#' Gradient descent whose step length is a scalar curvature estimate computed
#' from the previous step's secant pair rather than found by searching.
#'
#' @param criterion The stopping rule; see \code{\link{crit_any}}.
#' @param variant \code{"alternate"} (default), \code{"bb1"} or \code{"bb2"};
#'   see Details.
#' @param alpha0 The step length used on the first iteration, before there is a
#'   secant pair to estimate one from. Defaults to \code{1e-2}.
#' @param alpha_min,alpha_max Bounds on the step length. Defaults \code{1e-10}
#'   and \code{1e10}.
#' @param curv_tol The relative curvature threshold: a secant pair is refused
#'   when \eqn{s^\top y \le c \lVert s \rVert \lVert y \rVert}{s'y <= c |s| |y|}
#'   for \eqn{c} equal to \code{curv_tol}. Defaults to \code{1e-10}, which is
#'   the same relative test \code{\link{bfgs}} applies to the same quantity.
#' @param step Initial multiplier offered to the line search. Defaults to
#'   \code{1}, so the Barzilai-Borwein step is tried unaltered first.
#' @param line_search The acceptance test for a trial step. Defaults to
#'   \code{\link{nonmonotone}}; see Details.
#' @param maxit Maximum iterations. Defaults to 1000.
#' @param max_eval Maximum objective evaluations. Defaults to 20000.
#' @param verbose Report progress? Defaults to \code{FALSE}.
#' @param refresh Report every this many iterations. Defaults to 20.
#' @param keep_trace Store the iteration path? Defaults to \code{FALSE}.
#'
#' @details
#' Take the direction \eqn{-\alpha g} with
#' \deqn{\alpha_{BB1} = \frac{s^\top s}{s^\top y}, \qquad
#'       \alpha_{BB2} = \frac{s^\top y}{y^\top y},}
#' the two Rayleigh quotients of the secant pair \eqn{s = x_k - x_{k-1}},
#' \eqn{y = g_k - g_{k-1}}. Both estimate the inverse curvature along the
#' direction just travelled, so this is a quasi-Newton method that has discarded
#' everything except one scalar. On a quadratic, where the curvature is
#' constant, that scalar is exactly right.
#'
#' On a quadratic, where the curvature is constant, the estimate is exact and
#' the method converges in two iterations; on a general smooth objective it
#' typically needs more iterations than \code{\link{bfgs}} while storing a
#' single scalar instead of a matrix.
#'
#' \subsection{Variants}{
#' \code{"bb1"} and \code{"bb2"} are the two quotients above, and
#' \code{"alternate"}, the default, switches between them at each iteration:
#' they estimate the same curvature from opposite ends, and alternating them is
#' more robust than either alone.
#' }
#'
#' \subsection{Line search}{
#' The Barzilai-Borwein step is offered to the line search first and unaltered,
#' and backtracking occurs only when it fails the acceptance test. Because the
#' method makes progress through steps that may increase the objective
#' temporarily, the default acceptance test is \code{\link{nonmonotone}}, which
#' requires improvement over the maximum of the last \code{memory} values
#' rather than over the current one; a plain Armijo condition rejects exactly
#' the steps the method relies on and slows it considerably.
#' \code{\link{nonmonotone}} with \code{memory = 0} coincides with
#' \code{\link{armijo}}.
#' }
#'
#' \subsection{Refused secant pairs}{
#' A pair is used only if it reports positive curvature by a relative margin,
#' \eqn{s^\top y > c \lVert s \rVert \lVert y \rVert}{s'y > c |s| |y|} with
#' \eqn{c} the \code{curv_tol} argument -- the same test \code{\link{bfgs}}
#' applies. When a pair is refused, the step length is reset to
#' \eqn{1/\lVert g \rVert_\infty}{1/max|g|}, giving a trial displacement of
#' order one in the parameters. This reset depends on the current gradient
#' rather than on the step length being replaced or on a fixed constant, so it
#' can neither freeze the iteration at a too-short step nor produce a step the
#' backtracking cannot rescale. Steps outside \code{[alpha_min, alpha_max]} are
#' clamped, and both the reset and the clamp are recorded in the trace.
#' }
#'
#' @return An S7 object of class \code{Bb}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @references
#' Barzilai, J. and Borwein, J. M. (1988). Two-point step size gradient methods.
#' \emph{IMA Journal of Numerical Analysis} \strong{8}, 141--148.
#'
#' @examples
#' bb()
#'
#' f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
#' gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
#'                     200 * (p[2] - p[1]^2))
#' minimize(bb(), f, c(-1.2, 1), gr = gr)@par
#'
#' # two iterations on a quadratic: one secant pair determines the curvature
#' minimize(bb(), function(p) sum((p - c(1, 2))^2), c(0, 0),
#'          gr = function(p) 2 * (p - c(1, 2)))@iterations
#'
#' @seealso \code{\link{gd}}, \code{\link{cg}}, \code{\link{lbfgs}}
#' @export
bb <- function(criterion = crit_any(crit_grad(1e-8), crit_rel_obj(1e-12)),
               variant = c("alternate", "bb1", "bb2"),
               alpha0 = 1e-2, alpha_min = 1e-10, alpha_max = 1e10,
               curv_tol = 1e-10,
               step = 1, line_search = nonmonotone(),
               maxit = 1000, max_eval = 20000,
               verbose = FALSE, refresh = 20, keep_trace = FALSE) {
  variant <- match.arg(variant)
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  check_tol(alpha0)
  check_tol(alpha_min)
  check_tol(alpha_max)
  check_tol(curv_tol)
  if (alpha_min >= alpha_max) {
    stop("'alpha_min' must be strictly below 'alpha_max'.", call. = FALSE)
  }
  Bb(
    name = paste0("barzilai-borwein (", variant, ")"), criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    variant = variant, alpha0 = alpha0,
    alpha_min = alpha_min, alpha_max = alpha_max, curv_tol = curv_tol,
    step = step, line_search = line_search
  )
}


#' @title Minimise by Barzilai-Borwein
#' @name minimize.Bb
#' @description Runs \code{\link{bb}} on the objective.
#' @param optimizer A \code{Bb} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Bb) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "bb", variant = optimizer@variant,
                     alpha0 = optimizer@alpha0,
                     alpha_min = optimizer@alpha_min,
                     alpha_max = optimizer@alpha_max,
                     curv_tol = optimizer@curv_tol))
  }
