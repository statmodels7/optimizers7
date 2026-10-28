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
#' A `GradientDescent` object holds the initial step length and the line
#' search only, since steepest descent does not keep a model of the surface.
#' It is built by [gd()].
#'
#' @details
#' Beyond the seven properties of every optimizer, a `GradientDescent`
#' carries `step`, `line_search` and `stop_decrement`, shared with the other
#' line-search methods.
#'
#' @param step The initial step length offered to the line search.
#' @param line_search A [line_search()] object.
#' @param stop_decrement The largest decrease at which a stopping rule on a
#'   change may end the run as converged; see [gd()].
#'
#' @return An S7 object of class `GradientDescent` inheriting from
#'   [optimizer()], with `step`, `line_search` and `stop_decrement` beside the
#'   seven shared properties.
#'
#' @seealso [gd()] for the constructor, [Cg] and [Bb] for the two other
#'   first-order methods.
#' @name GradientDescent-class
#' @aliases GradientDescent
#' @keywords internal
GradientDescent <- S7::new_class("GradientDescent", parent = optimizer,
  properties = list(
    step        = S7::class_numeric,
    line_search = S7::class_any,
    stop_decrement = S7::class_numeric
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
#'   [armijo()] and [wolfe()]. Defaults to Armijo backtracking, which is
#'   sufficient for a method that keeps no curvature approximation.
#' @param maxit Maximum iterations. Defaults to 500.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`, which
#'   sets no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 10.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#' @param stop_decrement The largest decrease
#'   \eqn{\tfrac12\lvert g^\top s\rvert}, with \eqn{s} the last step and
#'   \eqn{g} the gradient at the new point, at which a stopping rule on the
#'   change in the parameters or in the objective, such as [crit_abs_par()] or
#'   [crit_abs_obj()], may end the run as converged. A single positive
#'   number, in the units of the objective, or `Inf` to accept every such
#'   stop. Defaults to `1e-3`. See below.
#'
#' @details
#' The direction \eqn{-g} minimizes \eqn{g^\top d} over directions of a given
#' Euclidean length. Under an exact line search consecutive directions are
#' orthogonal, so on an ill-conditioned objective the iterates zigzag and
#' converge slowly. [cg()] corrects this by combining each direction with
#' the previous one, at no extra cost per iteration.
#'
#' # A stop on a change
#'
#' A stopping rule on the change in the parameters or in the objective states
#' that the iterate no longer moves, and near a minimum that is the same as
#' being there. Away from one it is not: a step can be short because the
#' objective is steep, or because the direction came from a poor model of the
#' surface, and the decrease \eqn{\tfrac12\lvert g^\top s\rvert} that the
#' step still predicts at the new point is then above what a minimum leaves.
#' A stop on the change is therefore not accepted where that decrease exceeds
#' `stop_decrement`. The run ends as not converged, with `criterion_met`
#' reporting `"stopped without converging"` and the message giving the
#' decrease, and the trace records the safeguard `change stop refused`. A
#' rule on the gradient, such as [crit_grad()], is not affected. The check
#' uses quantities that the loop already holds. `stop_decrement = Inf`
#' accepts every stop on a change.
#'
#' Gradient descent serves as a baseline for the other methods. On the
#' Rosenbrock function from the usual start it does not converge within its
#' default budget of 500 iterations, while [cg()], [bb()] and [bfgs()]
#' converge in far fewer (the last example).
#'
#' @return An S7 object of class [GradientDescent], inheriting from
#'   [optimizer()], to be handed to [minimize()].
#'
#' @examples
#' gd()
#' gd(criterion = crit_grad(1e-10), maxit = 2000)
#'
#' # A well-conditioned problem is solved quickly.
#' minimize(gd(), function(p) sum((p - c(1, 2))^2), c(0, 0),
#'          gr = function(p) 2 * (p - c(1, 2)))
#'
#' # In a curved valley the budget runs out before convergence.
#' ros <- test_problems("rosenbrock")[[1]]
#' r <- minimize(gd(), ros$fn, ros$par, gr = ros$gr)
#' c(converged = r@converged, iterations = r@iterations,
#'   evaluations = r@counts[["f"]], value = r@value)
#'
#' @seealso [cg()] for the same storage and a better direction, [bb()] for a
#'   scalar curvature estimate, [bfgs()] for a matrix one.
#' @references
#' Cauchy's method is the oldest of the methods in the package; chapter 3 of
#' the following book treats it, including the dependence of its linear rate
#' on the condition number:
#' Nocedal, J. and Wright, S. J. (2006). *Numerical Optimization*,
#' 2nd edition. Springer, New York.
#'
#' @export
gd <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
               step = 1, line_search = armijo(),
               maxit = 500, max_eval = Inf,
               verbose = FALSE, refresh = 10, keep_trace = FALSE,
               stop_decrement = 1e-3) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  check_stop_decrement(stop_decrement)
  GradientDescent(
    name = "gradient descent", criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    step = step, line_search = line_search,
    stop_decrement = stop_decrement
  )
}


#' @title Minimize by Gradient Descent
#' @name minimize.GradientDescent
#'
#' @description
#' Runs [gd()] on the objective: the direction is \eqn{-g}, and the line
#' search chooses the step length. It shares the descent loop with [cg()]
#' and [bb()], which differ from it only in the direction.
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
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "gd", stop_decrement = optimizer@stop_decrement))
  }


# --- conjugate gradients ----------------------------------------------------

#' @title S7 Class for Conjugate Gradients
#'
#' @description
#' A `Cg` object holds which of the four \eqn{\beta} formulas combines the
#' gradient with the previous direction, how often the method is restarted
#' at steepest descent, and the step length and the line search. It is built
#' by [cg()].
#'
#' @details
#' Beyond the seven properties of every optimizer, a `Cg` carries five:
#' `step`, `line_search` and `stop_decrement`, shared with the other
#' line-search methods, and `beta` and `restart_every`, which are its own. The
#' previous direction is not a property; it exists only inside the run.
#'
#' @param beta Which update formula, one of `"pr"`, `"fr"`, `"hs"`, `"dy"`.
#' @param restart_every How often the method is restarted at steepest
#'   descent; `0` means never.
#' @param step The initial step length offered to the line search.
#' @param line_search A [line_search()] object.
#' @param stop_decrement The largest decrease at which a stopping rule on a
#'   change may end the run as converged; see [gd()].
#'
#' @return An S7 object of class `Cg` inheriting from [optimizer()], with the
#'   five properties above beside the seven shared ones.
#'
#' @seealso [cg()] for the constructor, [GradientDescent] for steepest
#'   descent.
#' @name Cg-class
#' @aliases Cg
#' @keywords internal
Cg <- S7::new_class("Cg", parent = optimizer,
  properties = list(
    beta          = S7::class_character,
    restart_every = S7::class_numeric,
    step          = S7::class_numeric,
    line_search   = S7::class_any,
    stop_decrement = S7::class_numeric
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
#' @param beta Which formula for \eqn{\beta}: `"pr"` (Polak-Ribière, the
#'   default), `"fr"` (Fletcher-Reeves), `"hs"` (Hestenes-Stiefel) or
#'   `"dy"` (Dai-Yuan). See Details.
#' @param restart_every Restart at steepest descent every this many iterations.
#'   Defaults to `0`, meaning never; a positive value is usually the
#'   dimension of the problem.
#' @param step Initial step length offered to the line search. Defaults to
#'   `1`.
#' @param line_search Defaults to `wolfe(c2 = 0.1)`; the convergence
#'   theory of the method assumes a Wolfe step, see Details.
#' @param maxit Maximum iterations. Defaults to 1000.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`, which
#'   sets no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 20.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#' @param stop_decrement The largest decrease
#'   \eqn{\tfrac12\lvert g^\top s\rvert} at which a stopping rule on the change
#'   in the parameters or in the objective may end the run as converged, as in
#'   [gd()]. A single positive number, in the units of the objective, or `Inf`
#'   to accept every such stop. Defaults to `1e-3`. See below.
#'
#' @details
#' The direction is
#' \deqn{d_k = -g_k + \beta_k d_{k-1},}
#' and the methods differ only in the choice of \eqn{\beta}. The storage is
#' two vectors against [bfgs()]'s \eqn{p \times p} matrix.
#'
#' # Conjugate directions and the line search
#'
#' On a quadratic with an exact line search the directions are conjugate
#' with respect to the Hessian, and the method terminates in \eqn{p} steps
#' without forming that Hessian. The line search here is inexact, so that
#' guarantee does not hold, and the number of iterations on a quadratic
#' depends on its conditioning. Each iteration still stores only two
#' vectors.
#'
#' # Choice of beta
#'
#' All four formulas agree on a quadratic under an exact line search and
#' differ elsewhere. `"fr"` has the simplest convergence theory and can
#' stall for many iterations after a poor step. `"pr"` recovers from a
#' poor step at once, since a small \eqn{y = g_k - g_{k-1}} makes \eqn{\beta}
#' small and the direction close to steepest descent; its known theoretical
#' non-convergence is removed by clamping \eqn{\beta} at zero, which restarts
#' the method and appears in the trace as `cg restart`. `"hs"` and `"dy"` are
#' the other two standard choices. `"pr"` is the default as the safest of the
#' four on a general objective; which formula needs the fewest iterations
#' depends on the problem (see the examples).
#'
#' # Line search
#'
#' The convergence results for these formulas assume a step that satisfies
#' the Wolfe conditions (the strong ones for `"fr"` and `"pr"`) and use them
#' to show that the direction produced is a descent direction. Backtracking
#' gives no such guarantee, so [wolfe()] is the default.
#'
#' [wolfe()] defaults to \eqn{c_2 = 0.9}, which suits [bfgs()], because its
#' curvature update compensates for an inexact step. The conjugacy
#' accumulated by conjugate gradients depends directly on the accuracy of the
#' line search, so the default here is `wolfe(c2 = 0.1)`.
#'
#' A direction that is not a descent direction is replaced by \eqn{-g}, and
#' the trace records the replacement as `cg not descent`. A periodic restart
#' set by `restart_every` is recorded as `cg restart`.
#'
#' # A stop on a change
#'
#' A stop on the change in the parameters or in the objective is not accepted
#' where the decrease \eqn{\tfrac12\lvert g^\top s\rvert} that the last step
#' still predicts at the new point exceeds `stop_decrement`; the run then ends
#' as not converged, as described in [gd()]. `stop_decrement = Inf` accepts
#' every such stop.
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
#' # The four formulas are compared on the same problem. They agree on a quadratic under
#' # an exact line search and differ here; which is best depends on the
#' # problem.
#' vapply(c("pr", "fr", "hs", "dy"),
#'        function(b) minimize(cg(beta = b), f, c(-1.2, 1), gr = gr)@iterations,
#'        integer(1))
#'
#' # Polak-Ribiere clamps beta at zero after a poor step, which restarts the
#' # method; the trace counts the restarts.
#' r <- minimize(cg(keep_trace = TRUE), f, c(-1.2, 1), gr = gr)
#' table(r@trace$safeguard)
#'
#' @seealso [gd()] for steepest descent, [lbfgs()] for a method with storage
#'   of the same order and more curvature information, [bb()] for the scalar
#'   estimate.
#' @export
cg <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
               beta = c("pr", "fr", "hs", "dy"), restart_every = 0,
               step = 1, line_search = wolfe(c2 = 0.1),
               maxit = 1000, max_eval = Inf,
               verbose = FALSE, refresh = 20, keep_trace = FALSE,
               stop_decrement = 1e-3) {
  beta <- match.arg(beta)
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  check_stop_decrement(stop_decrement)
  if (length(restart_every) != 1L || !is.numeric(restart_every) ||
      is.na(restart_every) || restart_every < 0) {
    stop("'restart_every' must be a single non-negative number.", call. = FALSE)
  }
  Cg(
    name = paste0("conjugate gradients (", beta, ")"), criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    beta = beta, restart_every = restart_every,
    step = step, line_search = line_search,
    stop_decrement = stop_decrement
  )
}


#' @title Minimize by Conjugate Gradients
#' @name minimize.Cg
#'
#' @description
#' Runs [cg()] on the objective: each direction is combined with the
#' previous one through the chosen \eqn{\beta}, and the line search chooses
#' the step length. It shares the descent loop with [gd()] and [bb()].
#'
#' @param optimizer A `Cg` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `he` is accepted
#'   and ignored. Bounds are taken and removed by reparametrization.
#'
#' @return An [optimizer_result()] whose trace, when kept, reports
#'   `cg restart` at each iteration where \eqn{\beta} was set to zero, because
#'   the Polak-Ribière value was negative or because `restart_every`
#'   iterations had passed, and `cg not descent` at each iteration where the
#'   direction was replaced by the negative gradient.
#'
#' @keywords internal
S7::method(minimize, Cg) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    run_descent(optimizer, fn, par, gr, he, lower, upper,
                list(type = "cg", beta = optimizer@beta,
                     restart_every = as.integer(optimizer@restart_every),
                     stop_decrement = optimizer@stop_decrement))
  }


# --- Barzilai-Borwein -------------------------------------------------------

#' @title S7 Class for the Barzilai-Borwein Method
#'
#' @description
#' A `Bb` object holds which of the two Rayleigh quotients gives the step
#' length, the bounds to which that quotient is clamped, the curvature
#' threshold that a secant pair must meet, and the step multiplier and the
#' line search. It is built by [bb()]. Its line search defaults to
#' [nonmonotone()]; it is the only shipped method with that default.
#'
#' @details
#' Beyond the seven properties of every optimizer, a `Bb` carries eight of
#' its own: `variant`, `alpha0`, `alpha_min`, `alpha_max` and `curv_tol` for
#' the step-length estimate, `step` and `line_search` for the search along the
#' direction, and `stop_decrement`.
#'
#' @param variant Which step-length formula, one of `"alternate"`, `"bb1"`,
#'   `"bb2"`.
#' @param alpha0 The step length used before there is a secant pair.
#' @param alpha_min,alpha_max The bounds to which the step length is
#'   clamped.
#' @param curv_tol The relative threshold below which a secant pair is
#'   rejected.
#' @param step,line_search The multiplier offered to the line search and the
#'   [line_search()] object, as in [bb()].
#' @param stop_decrement The largest decrease at which a stopping rule on a
#'   change may end the run as converged; see [gd()].
#'
#' @return An S7 object of class `Bb` inheriting from [optimizer()], with the
#'   eight properties above beside the seven shared ones.
#'
#' @seealso [bb()] for the constructor, [nonmonotone()] for the acceptance
#'   test that is its default.
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
    line_search = S7::class_any,
    stop_decrement = S7::class_numeric
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
#'   secant pair from which to estimate one. Defaults to `1e-2`.
#' @param alpha_min,alpha_max Bounds on the step length. The defaults are
#'   `1e-10` and `1e10`.
#' @param curv_tol The relative curvature threshold: a secant pair is rejected
#'   when \eqn{s^\top y \le c \lVert s \rVert \lVert y \rVert}{s'y <= c |s| |y|}
#'   for \eqn{c} equal to `curv_tol`. Defaults to `1e-10`, which is
#'   the same relative test that [bfgs()] applies to the same quantity.
#' @param step Initial multiplier offered to the line search. Defaults to
#'   `1`, so the Barzilai-Borwein step is tried unaltered first.
#' @param line_search The acceptance test for a trial step. Defaults to
#'   [nonmonotone()]; see Details.
#' @param maxit Maximum iterations. Defaults to 1000.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`, which
#'   sets no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 20.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#' @param stop_decrement The largest decrease
#'   \eqn{\tfrac12\lvert g^\top s\rvert} at which a stopping rule on the change
#'   in the parameters or in the objective may end the run as converged, as in
#'   [gd()]. A single positive number, in the units of the objective, or `Inf`
#'   to accept every such stop. Defaults to `1e-3`. See below.
#'
#' @details
#' Take the direction \eqn{-\alpha g} with
#' \deqn{\alpha_{BB1} = \frac{s^\top s}{s^\top y}, \qquad
#'       \alpha_{BB2} = \frac{s^\top y}{y^\top y},}
#' the two Rayleigh quotients of the secant pair \eqn{s = x_k - x_{k-1}},
#' \eqn{y = g_k - g_{k-1}}. Both estimate the inverse curvature along the
#' direction just traveled, so the method is a quasi-Newton method that keeps
#' a single scalar. When the Hessian is a multiple of the identity, the
#' curvature is the same in every direction, the estimate from one secant
#' pair is exact, and the method converges in two iterations. Otherwise
#' the scalar is the curvature along the last step only, and the number of
#' iterations grows with the spread of the curvatures. On a general smooth
#' objective the method usually needs more iterations and more evaluations
#' than [bfgs()], and it stores one scalar instead of a matrix.
#'
#' # Variants
#'
#' `"bb1"` and `"bb2"` are the two quotients above, and `"alternate"`, the
#' default, switches between them at each iteration. They estimate the same
#' curvature in two ways, and alternating is the more robust choice across
#' problems, although not always the fastest on a given one (see the
#' examples).
#'
#' # Line search
#'
#' The Barzilai-Borwein step is offered to the line search first and
#' unaltered, and backtracking occurs only when it fails the acceptance test.
#' The method makes progress through steps that may increase the objective
#' temporarily, so the default acceptance test is [nonmonotone()], which
#' requires an improvement over the largest of the last `memory + 1`
#' objective values, which include the current one, and not only over the
#' current value. A plain Armijo condition rejects every step that increases the
#' objective, and with it the steps on which the method relies, so the run
#' needs more than twice as many evaluations (see the examples).
#'
#' `nonmonotone(memory = 0)` is [armijo()] value for value, so a comparison
#' of the two isolates the effect of the memory.
#'
#' # Rejected secant pairs
#'
#' A pair is used only if it reports positive curvature by a relative margin,
#' \eqn{s^\top y > c \lVert s \rVert \lVert y \rVert}{s'y > c |s| |y|} with
#' \eqn{c} the `curv_tol` argument, the same test that [bfgs()] applies. When a
#' pair is rejected the step length is reset to
#' \eqn{1/\lVert g \rVert_\infty}{1/max|g|}, a trial displacement of order one
#' in the parameters. The reset reads the current gradient and not the step
#' length being replaced or a fixed constant, so it cannot freeze the
#' iteration at a too-short step and cannot produce a step that the
#' backtracking is unable to rescale. Steps outside `[alpha_min, alpha_max]`
#' are clamped. The reset and the clamp appear in the trace as
#' `bb curvature reset` and `bb step clamped`, and a step shortened by the
#' line search as `step shortened`.
#'
#' # A stop on a change
#'
#' A stop on the change in the parameters or in the objective is not accepted
#' where the decrease \eqn{\tfrac12\lvert g^\top s\rvert} that the last step
#' still predicts at the new point exceeds `stop_decrement`; the run then ends
#' as not converged, as described in [gd()]. `stop_decrement = Inf` accepts
#' every such stop.
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
#' # Two iterations on a quadratic whose Hessian is a multiple of the
#' # identity: one secant pair determines the curvature exactly.
#' minimize(bb(), function(p) sum((p - c(1, 2))^2), c(0, 0),
#'          gr = function(p) 2 * (p - c(1, 2)))@iterations
#'
#' # The nonmonotone search against the monotone one. With memory = 0 it is
#' # armijo() value for value, and both give the same, slower run.
#' evals <- function(ls) minimize(bb(line_search = ls), f, c(-1.2, 1),
#'                                gr = gr)@counts[["f"]]
#' c(nonmonotone = evals(nonmonotone()),
#'   memory_zero = evals(nonmonotone(memory = 0)),
#'   armijo      = evals(armijo()))
#'
#' # Which variant is fastest depends on the problem.
#' vapply(c("bb1", "bb2", "alternate"),
#'        function(v) minimize(bb(variant = v), f, c(-1.2, 1), gr = gr)@iterations,
#'        integer(1))
#'
#' @seealso [gd()] for the same direction with a line-searched step,
#'   [nonmonotone()] for the default acceptance test of this method,
#'   [lbfgs()] for a method that keeps more curvature information.
#' @export
bb <- function(criterion = crit_any(crit_grad(), crit_abs_obj(), crit_abs_par()),
               variant = c("alternate", "bb1", "bb2"),
               alpha0 = 1e-2, alpha_min = 1e-10, alpha_max = 1e10,
               curv_tol = 1e-10,
               step = 1, line_search = nonmonotone(),
               maxit = 1000, max_eval = Inf,
               verbose = FALSE, refresh = 20, keep_trace = FALSE,
               stop_decrement = 1e-3) {
  variant <- match.arg(variant)
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_step(step)
  check_line_search(line_search)
  check_stop_decrement(stop_decrement)
  check_tol(alpha0, "alpha0")
  check_tol(alpha_min, "alpha_min")
  check_tol(alpha_max, "alpha_max")
  check_tol(curv_tol, "curv_tol")
  if (alpha_min >= alpha_max) {
    stop("'alpha_min' must be strictly below 'alpha_max'.", call. = FALSE)
  }
  Bb(
    name = paste0("barzilai-borwein (", variant, ")"), criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    variant = variant, alpha0 = alpha0,
    alpha_min = alpha_min, alpha_max = alpha_max, curv_tol = curv_tol,
    step = step, line_search = line_search,
    stop_decrement = stop_decrement
  )
}


#' @title Minimize by Barzilai-Borwein
#' @name minimize.Bb
#'
#' @description
#' Runs [bb()] on the objective: the step length is formed from the
#' previous secant pair and offered to the line search unaltered, which
#' backtracks only if it is rejected. It shares the descent loop with [gd()]
#' and [cg()].
#'
#' @param optimizer A `Bb` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `he` is accepted
#'   and ignored. Bounds are taken and removed by reparametrization.
#'
#' @return An [optimizer_result()] whose trace, when kept, reports
#'   `bb curvature reset` where a secant pair had no usable curvature,
#'   `bb step clamped` where the step was clamped to
#'   `[alpha_min, alpha_max]`, and `step shortened` where the line search
#'   backtracked.
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
                     curv_tol = optimizer@curv_tol,
                     stop_decrement = optimizer@stop_decrement))
  }
