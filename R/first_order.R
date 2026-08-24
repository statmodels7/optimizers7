#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

# The methods that use the gradient and nothing else. All three are a direction
# plus a line search, so they share the compiled loop of descent.cpp and differ
# only in how the direction is formed.


#' @title S7 Class for Gradient Descent
#'
#' @description
#' An optimizer holding the initial step length and the line search, and
#' nothing else: steepest descent carries no model of the surface, so there
#' is nothing else to hold. Built by [gd()].
#'
#' @details
#' Beyond the seven properties every optimizer has, a `GradientDescent`
#' carries `step` and `line_search`, both shared with the other line-search
#' methods. It is the smallest optimizer class in the package.
#'
#' @param step The initial step length offered to the line search.
#' @param line_search A [line_search()] object.
#'
#' @return An S7 object of class `GradientDescent` inheriting from
#'   [optimizer()], with `step` and `line_search` beside the seven shared
#'   properties.
#'
#' @seealso [gd()] for the constructor, [Cg] and [Bb] for the two first-order
#'   methods that carry a little more.
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
#' Steepest descent: the search direction is the negative gradient
#' \eqn{-g(x)}, and a line search chooses the step length.
#'
#' @param criterion The stopping rule; see [crit_any()].
#' @param step Initial step length offered to the line search each iteration.
#'   Defaults to `1`.
#' @param line_search How far to go along the direction; see
#'   [armijo()] and [wolfe()]. Defaults to Armijo
#'   backtracking, which suffices for a method carrying no curvature
#'   approximation to protect.
#' @param maxit Maximum iterations. Defaults to 500.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`:
#'   no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 10.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#'
#' @details
#' The direction \eqn{-g} minimizes \eqn{g^\top d} over directions of a given
#' Euclidean length. Under an exact line search consecutive directions are
#' orthogonal, so on an ill-conditioned objective the iterates zigzag and
#' converge slowly. [cg()] repairs that by combining each direction with the
#' previous one at no extra cost per iteration.
#'
#' Gradient descent is a baseline. On Rosenbrock from the customary start it
#' exhausts its budget of 500 iterations at a value of `1.9e-03` after 4905
#' objective evaluations, where [cg()] reaches `7.4e-10` in 25 iterations,
#' [bb()] `5.2e-09` in 58 and [bfgs()] `3.2e-13` in 35.
#'
#' @return An S7 object of class [GradientDescent], inheriting from
#'   [optimizer()], to be handed to [minimize()].
#'
#' @examples
#' gd()
#' gd(criterion = crit_grad(1e-10), maxit = 2000)
#'
#' # It solves an easy problem readily.
#' minimize(gd(), function(p) sum((p - c(1, 2))^2), c(0, 0),
#'          gr = function(p) 2 * (p - c(1, 2)))
#'
#' # And is the wrong tool for a curved valley: the budget runs out first.
#' ros <- test_problems("rosenbrock")[[1]]
#' r <- minimize(gd(), ros$fn, ros$par, gr = ros$gr)
#' c(converged = r@converged, iterations = r@iterations,
#'   evaluations = r@counts[["f"]], value = r@value)
#'
#' @seealso [cg()] for the same storage and a better direction, [bb()] for a
#'   scalar curvature estimate, [bfgs()] for a matrix one.
#' @references
#' Cauchy's method is the oldest of them; the modern treatment,
#' including why its rate is linear in the condition number, is
#' chapter 3 of
#' Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*,
#' 2nd edition. Springer, New York.
#'
#' @export
gd <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
               step = 1, line_search = armijo(),
               maxit = 500, max_eval = Inf,
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


#' @title Minimize by Gradient Descent
#' @name minimize.GradientDescent
#'
#' @description
#' Runs [gd()] on the objective: take \eqn{-g} as the direction and let the
#' line search choose how far. Shares the descent loop with [cg()] and
#' [bb()], which differ from it only in how the direction is formed.
#'
#' @param optimizer A `GradientDescent` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `he` is accepted
#'   and ignored. Bounds are taken and removed by reparametrization.
#'
#' @return An [optimizer_result()] whose `gradient` is \eqn{\nabla f} at
#'   `par` and whose trace, when kept, carries a `gnorm` column.
#'
#' @keywords internal
S7::method(minimize, GradientDescent) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    run_descent(optimizer, fn, par, gr, he, lower, upper, list(type = "gd"))
  }


# --- conjugate gradients ----------------------------------------------------

#' @title S7 Class for Conjugate Gradients
#'
#' @description
#' An optimizer holding which of the four \eqn{\beta} formulas bends the
#' direction, how often the method is restarted at steepest descent, and the
#' usual step length and line search. Built by [cg()].
#'
#' @details
#' Beyond the seven properties every optimizer has, a `Cg` carries four:
#' `step` and `line_search`, shared with the other line-search methods, and
#' `beta` and `restart_every`, which are its own. The previous direction is
#' not a property; it lives in the run.
#'
#' @param beta Which update formula, one of `"pr"`, `"fr"`, `"hs"`, `"dy"`.
#' @param restart_every How often the method is restarted at steepest
#'   descent; `0` means never.
#' @param step The initial step length offered to the line search.
#' @param line_search A [line_search()] object.
#'
#' @return An S7 object of class `Cg` inheriting from [optimizer()], with the
#'   four properties above beside the seven shared ones.
#'
#' @seealso [cg()] for the constructor, [GradientDescent] for the method it
#'   improves on.
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
#' Nonlinear conjugate gradients: the direction is
#' \eqn{d_k = -g_k + \beta_k d_{k-1}}, so consecutive directions are conjugate
#' on a quadratic rather than orthogonal. Storage is two vectors, with no
#' matrix, which suits problems too large for a Hessian.
#'
#' @param criterion The stopping rule; see [crit_any()].
#' @param beta Which formula for the bend: `"pr"` (Polak-Ribière, the
#'   default), `"fr"` (Fletcher-Reeves), `"hs"` (Hestenes-Stiefel) or
#'   `"dy"` (Dai-Yuan). See Details.
#' @param restart_every Restart at steepest descent every this many iterations.
#'   Defaults to `0`, meaning never; a positive value is usually the
#'   dimension of the problem.
#' @param step Initial step length offered to the line search. Defaults to
#'   `1`.
#' @param line_search Defaults to `wolfe(c2 = 0.1)`; the convergence
#'   theory of the method assumes a strong Wolfe step, see Details.
#' @param maxit Maximum iterations. Defaults to 1000.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`:
#'   no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 20.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#'
#' @details
#' The direction is
#' \deqn{d_k = -g_k + \beta_k d_{k-1},}
#' and the methods differ only in the choice of \eqn{\beta}. The storage is
#' two vectors against [bfgs()]'s \eqn{p \times p} matrix.
#'
#' # What conjugacy buys, and what it needs
#'
#' On a quadratic **with an exact line search** the directions come out
#' conjugate with respect to the Hessian, and the method then terminates in
#' \eqn{p} steps without ever forming that Hessian. The line search here is
#' inexact, so that guarantee does not transfer: measured on dense quadratics
#' at `crit_grad(1e-12)`, the run takes 16, 34 and 70 iterations at
#' \eqn{p = 3, 8, 20}. The saving is still real, since each iteration costs
#' two vectors.
#'
#' # Choice of beta
#'
#' All four agree on a quadratic under an exact line search and differ
#' everywhere else. `"fr"` has the cleanest convergence theory and the
#' well-known practical fault of stalling for many iterations after a poor
#' step. `"pr"` recovers from a poor step immediately, a small
#' \eqn{y = g_k - g_{k-1}} sending \eqn{\beta} towards zero and the method
#' back to steepest descent; its known theoretical non-convergence is
#' repaired by clamping \eqn{\beta} at zero, which restarts the method and
#' appears in the trace as `cg restart`. `"hs"` and `"dy"` are the other two
#' standard choices.
#'
#' Measured on Rosenbrock from the customary start, iterations and objective
#' evaluations: `pr` 25 and 249 with twelve clamps, `fr` 59 and 743, `hs` 17
#' and 166, `dy` 50 and 667. `pr` is the default as the safest of the four on
#' a general objective, and `hs` was the fastest here; on another problem the
#' order will differ.
#'
#' # Line search
#'
#' The theory behind every one of these formulas assumes a step satisfying
#' the **strong** Wolfe conditions, and uses it to prove that the direction
#' produced is a descent direction at all. Backtracking gives no such
#' guarantee, so [wolfe()] is the default and departing from it departs from
#' the theory.
#'
#' The constant matters too. [wolfe()] defaults to \eqn{c_2 = 0.9}, which is
#' right for [bfgs()], where the curvature approximation repairs a loose
#' step. Conjugate gradients has nothing to repair with: the accumulated
#' conjugacy is only as good as the line search that produced it. The default
#' here is therefore `wolfe(c2 = 0.1)`, and the cost of loosening it is
#' modest on Rosenbrock, 28 iterations against 25, with the value reached
#' `3.3e-08` against `7.4e-10`.
#'
#' A direction that comes out non-descent is replaced by \eqn{-g} and the
#' substitution is reported in the trace, as a safeguard against the cases
#' the theory misses.
#'
#' @return An S7 object of class [Cg], inheriting from [optimizer()], to be
#'   handed to [minimize()].
#'
#' @references
#' Hestenes, M. R. and Stiefel, E. (1952). Methods of conjugate gradients for
#' solving linear systems. *Journal of Research of the NBS* **49**,
#' 409--436.
#'
#' Polak, E. and Ribière, G. (1969). Note sur la convergence de méthodes de
#' directions conjuguées. *Revue Française d'Informatique et de Recherche
#' Opérationnelle* **3**, 35--43.
#'
#' Dai, Y. H. and Yuan, Y. (1999). A nonlinear conjugate gradient method with a
#' strong global convergence property. *SIAM Journal on Optimization*
#' **10**, 177--182.
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
#' # The four formulas on the same problem. They agree on a quadratic under an
#' # exact line search and differ here, and which is best is a property of the
#' # problem.
#' vapply(c("pr", "fr", "hs", "dy"),
#'        function(b) minimize(cg(beta = b), f, c(-1.2, 1), gr = gr)@iterations,
#'        integer(1))
#'
#' # Polak-Ribiere clamps beta at zero when a step goes badly, which restarts
#' # the method; the trace counts those.
#' r <- minimize(cg(keep_trace = TRUE), f, c(-1.2, 1), gr = gr)
#' table(r@trace$safeguard)
#'
#' @seealso [gd()] for the direction this one bends, [lbfgs()] for a method
#'   with the same storage order and more curvature, [bb()] for the scalar
#'   estimate.
#' @export
cg <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
               beta = c("pr", "fr", "hs", "dy"), restart_every = 0,
               step = 1, line_search = wolfe(c2 = 0.1),
               maxit = 1000, max_eval = Inf,
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


#' @title Minimize by Conjugate Gradients
#' @name minimize.Cg
#'
#' @description
#' Runs [cg()] on the objective: bend each direction with the previous one by
#' the chosen \eqn{\beta}, and let the line search choose how far along it to
#' go. Shares the descent loop with [gd()] and [bb()].
#'
#' @param optimizer A `Cg` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `he` is accepted
#'   and ignored. Bounds are taken and removed by reparametrization.
#'
#' @return An [optimizer_result()] whose trace, when kept, reports
#'   `cg restart` at each iteration where the bend was clamped to zero or the
#'   direction was replaced by the gradient.
#'
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
#'
#' @description
#' An optimizer holding which of the two Rayleigh quotients gives the step
#' length, the bounds that quotient is clamped to, the curvature threshold a
#' secant pair must meet, and the usual step multiplier and line search.
#' Built by [bb()]. Its line search defaults to [nonmonotone()], the only
#' shipped method for which it does.
#'
#' @details
#' Beyond the seven properties every optimizer has, a `Bb` carries seven of
#' its own: `variant`, `alpha0`, `alpha_min`, `alpha_max` and `curv_tol` for
#' the step-length estimate, and `step` and `line_search` for what is done
#' with it.
#'
#' @param variant Which step-length formula, one of `"alternate"`, `"bb1"`,
#'   `"bb2"`.
#' @param alpha0 The step length used before there is a secant pair.
#' @param alpha_min,alpha_max Bounds the step length is clamped to.
#' @param curv_tol The relative threshold below which a secant pair is
#'   rejected.
#' @param step,line_search The multiplier offered to the line search and the
#'   [line_search()] object, as in [bb()].
#'
#' @return An S7 object of class `Bb` inheriting from [optimizer()], with the
#'   seven properties above beside the seven shared ones.
#'
#' @seealso [bb()] for the constructor, [nonmonotone()] for the acceptance
#'   test it defaults to.
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
#' Gradient descent with the Barzilai-Borwein step length: the direction is
#' \eqn{-\alpha_k g_k}, where \eqn{\alpha_k} is a scalar estimate of the
#' inverse curvature computed from the previous secant pair.
#'
#' @param criterion The stopping rule; see [crit_any()].
#' @param variant `"alternate"` (default), `"bb1"` or `"bb2"`;
#'   see Details.
#' @param alpha0 The step length used on the first iteration, before there is a
#'   secant pair to estimate one from. Defaults to `1e-2`.
#' @param alpha_min,alpha_max Bounds on the step length. Defaults `1e-10`
#'   and `1e10`.
#' @param curv_tol The relative curvature threshold: a secant pair is rejected
#'   when \eqn{s^\top y \le c \lVert s \rVert \lVert y \rVert}{s'y <= c |s| |y|}
#'   for \eqn{c} equal to `curv_tol`. Defaults to `1e-10`, which is
#'   the same relative test [bfgs()] applies to the same quantity.
#' @param step Initial multiplier offered to the line search. Defaults to
#'   `1`, so the Barzilai-Borwein step is tried unaltered first.
#' @param line_search The acceptance test for a trial step. Defaults to
#'   [nonmonotone()]; see Details.
#' @param maxit Maximum iterations. Defaults to 1000.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`:
#'   no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 20.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#'
#' @details
#' Take the direction \eqn{-\alpha g} with
#' \deqn{\alpha_{BB1} = \frac{s^\top s}{s^\top y}, \qquad
#'       \alpha_{BB2} = \frac{s^\top y}{y^\top y},}
#' the two Rayleigh quotients of the secant pair \eqn{s = x_k - x_{k-1}},
#' \eqn{y = g_k - g_{k-1}}. Both estimate the inverse curvature along the
#' direction just traveled, so this is a quasi-Newton method that has discarded
#' everything except one scalar. On a quadratic, where the curvature is
#' constant, that scalar is exactly right.
#'
#' On a quadratic, where the curvature is constant, the estimate is exact and
#' the method converges in **two** iterations. On a general smooth objective
#' it needs more iterations than [bfgs()] while storing a single scalar
#' instead of a matrix: on Rosenbrock 58 iterations against 35, but 67
#' objective evaluations against 49.
#'
#' # Variants
#'
#' `"bb1"` and `"bb2"` are the two quotients above, and `"alternate"`, the
#' default, switches between them at each iteration. They estimate the same
#' curvature from opposite ends, and alternating is the more robust choice
#' across problems, though not always the fastest on any one: measured on
#' Rosenbrock, `bb1` takes 56 iterations and 97 evaluations, `bb2` 51 and 56,
#' `alternate` 58 and 67.
#'
#' # Line search
#'
#' The Barzilai-Borwein step is offered to the line search first and
#' unaltered, and backtracking occurs only when it fails the acceptance test.
#' The method makes progress through steps that may increase the objective
#' temporarily, so the default acceptance test is [nonmonotone()], which asks
#' for improvement over the maximum of the last `memory` values instead of
#' over the current one. A plain Armijo condition rejects exactly the steps
#' the method relies on: measured on Rosenbrock, 72 iterations and 154
#' evaluations against 58 and 67.
#'
#' `nonmonotone(memory = 0)` is [armijo()] value for value, and the two give
#' the identical run here, 72 iterations and 154 evaluations. That identity is
#' what makes the comparison a comparison of the memory alone.
#'
#' # Rejected secant pairs
#'
#' A pair is used only if it reports positive curvature by a relative margin,
#' \eqn{s^\top y > c \lVert s \rVert \lVert y \rVert}{s'y > c |s| |y|} with
#' \eqn{c} the `curv_tol` argument, the same test [bfgs()] applies. When a
#' pair is rejected the step length is reset to
#' \eqn{1/\lVert g \rVert_\infty}{1/max|g|}, a trial displacement of order one
#' in the parameters. The reset reads the current gradient and not the step
#' length being replaced or a fixed constant, so it can neither freeze the
#' iteration at a too-short step nor produce one the backtracking cannot
#' rescale. Steps outside `[alpha_min, alpha_max]` are clamped, and both the
#' reset and the clamp appear in the trace, as `bb curvature reset` and
#' `step shortened`.
#'
#' @return An S7 object of class [Bb], inheriting from [optimizer()], to be
#'   handed to [minimize()].
#'
#' @references
#' Barzilai, J. and Borwein, J. M. (1988). Two-point step size gradient methods.
#' *IMA Journal of Numerical Analysis* **8**, 141--148.
#'
#' @examples
#' bb()
#'
#' f  <- function(p) 100 * (p[2] - p[1]^2)^2 + (1 - p[1])^2
#' gr <- function(p) c(-400 * p[1] * (p[2] - p[1]^2) - 2 * (1 - p[1]),
#'                     200 * (p[2] - p[1]^2))
#' minimize(bb(), f, c(-1.2, 1), gr = gr)@par
#'
#' # Two iterations on a quadratic: one secant pair determines the curvature,
#' # and on a quadratic that curvature is exactly right.
#' minimize(bb(), function(p) sum((p - c(1, 2))^2), c(0, 0),
#'          gr = function(p) 2 * (p - c(1, 2)))@iterations
#'
#' # The nonmonotone rule is what the method needs. With memory = 0 it is
#' # armijo() value for value, and both give the identical, slower run.
#' evals <- function(ls) minimize(bb(line_search = ls), f, c(-1.2, 1),
#'                                gr = gr)@counts[["f"]]
#' c(nonmonotone = evals(nonmonotone()),
#'   memory_zero = evals(nonmonotone(memory = 0)),
#'   armijo      = evals(armijo()))
#'
#' # Which variant is fastest is a property of the problem.
#' vapply(c("bb1", "bb2", "alternate"),
#'        function(v) minimize(bb(variant = v), f, c(-1.2, 1), gr = gr)@iterations,
#'        integer(1))
#'
#' @seealso [gd()] for the same direction with a line-searched step,
#'   [nonmonotone()] for the acceptance test this method needs, [lbfgs()]
#'   for the next amount of curvature to carry.
#' @export
bb <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
               variant = c("alternate", "bb1", "bb2"),
               alpha0 = 1e-2, alpha_min = 1e-10, alpha_max = 1e10,
               curv_tol = 1e-10,
               step = 1, line_search = nonmonotone(),
               maxit = 1000, max_eval = Inf,
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


#' @title Minimize by Barzilai-Borwein
#' @name minimize.Bb
#'
#' @description
#' Runs [bb()] on the objective: form the step length from the previous
#' secant pair, offer that step to the line search unaltered, and backtrack
#' only if it is rejected. Shares the descent loop with [gd()] and [cg()].
#'
#' @param optimizer A `Bb` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `he` is accepted
#'   and ignored. Bounds are taken and removed by reparametrization.
#'
#' @return An [optimizer_result()] whose trace, when kept, reports
#'   `bb curvature reset` where a secant pair carried none and
#'   `step shortened` where the step was clamped or backtracked.
#'
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
