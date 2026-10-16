#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

#' @title S7 Class for the Proximal Bundle Method
#'
#' @description
#' This class is an optimizer that holds the proximal parameter \eqn{t} and
#' its bounds, the acceptance fraction that separates a serious step from a
#' null one, the size of the bundle of linearizations, and the effort spent on
#' the subproblem. It is built by [bundle()].
#'
#' @details
#' Beyond the seven properties of every optimizer, a `Bundle` carries seven
#' of its own, in three groups: the proximal parameter (`t0`, `t_min`,
#' `t_max`), the acceptance test (`m_serious`), and the model and its
#' subproblem (`bundle_size`, `qp_iters`, `qp_tol`).
#'
#' @param t0 Initial value of the proximal parameter \eqn{t}, as a step
#'   length on the parameter scale.
#' @param t_min,t_max Bounds on \eqn{t}.
#' @param m_serious Fraction of the predicted decrease that a serious step
#'   must achieve.
#' @param bundle_size Largest number of linearizations kept.
#' @param qp_iters,qp_tol Effort spent on the subproblem.
#'
#' @return An S7 object of class `Bundle` inheriting from [optimizer()], with
#'   the seven properties above beside the seven shared ones.
#'
#' @seealso [bundle()] for the constructor, [crit_stationary()] for its
#'   default rule.
#' @name Bundle-class
#' @aliases Bundle
#' @keywords internal
Bundle <- S7::new_class("Bundle", parent = optimizer,
  properties = list(
    t0          = S7::class_numeric,
    t_min       = S7::class_numeric,
    t_max       = S7::class_numeric,
    m_serious   = S7::class_numeric,
    bundle_size = S7::class_numeric,
    qp_iters    = S7::class_numeric,
    qp_tol      = S7::class_numeric
  ))


#' @title The Proximal Bundle Method
#'
#' @description
#' Proximal bundle method for convex non-smooth objectives: subgradients
#' collected at the visited points define a cutting-plane model, a proximal
#' term keeps the step near the stability center, and convergence is tested
#' on the optimality estimate \eqn{\lVert p \rVert^2 + \alpha}, built from
#' the aggregate subgradient \eqn{p} and its linearization error
#' \eqn{\alpha}, and not on a single subgradient, which does not vanish at a
#' kink.
#'
#' @param criterion The stopping rule. Defaults to `crit_stationary(1e-8)`,
#'   read on the optimality estimate \eqn{\lVert p \rVert^2 + \alpha}.
#' @param t0 Initial value of the proximal parameter \eqn{t}, expressed as a
#'   step length on the parameter scale and not as a bare multiplier;
#'   see Details. Defaults to `1`.
#' @param t_min,t_max Bounds on it, so that a run of null steps or a run of
#'   serious ones cannot drive it to zero or to infinity. The defaults are
#'   `1e-10` and `1e10`.
#' @param m_serious Fraction of the predicted decrease that a step must actually
#'   deliver to be accepted, in \eqn{(0, 1)}. Defaults to `0.1`.
#' @param bundle_size Largest number of linearizations kept; when a new one
#'   would exceed it, the oldest are replaced by the aggregate. Defaults to
#'   `20`.
#' @param qp_iters,qp_tol Effort spent on the subproblem: at most this many
#'   accelerated projected-gradient steps, stopping when the weights move by
#'   less than `qp_tol`. The defaults are `500` and `1e-12`.
#' @param maxit Maximum iterations. Defaults to 500.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`, which
#'   sets no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 10.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#'
#' @details
#' # Behavior at a kink
#'
#' Near the minimum of \eqn{\lvert x \rvert} the evaluated subgradient is
#' \eqn{\pm 1}. The gradient norm seen by a descent method therefore stays
#' large, the line search finds no acceptable step, and the run stops at the
#' minimizer without meeting its stopping rule, whatever the tolerance. On
#' the median example below, [bfgs()] reaches the median and reports
#' `converged = FALSE`.
#'
#' The bundle method tests no single subgradient. It keeps a *collection* of
#' them, from the points that it has visited, and builds the piecewise-linear
#' model \eqn{\max_j \{ f_j + g_j'(x - x_j) \}}, the largest of the
#' linearizations. Such a model represents a kink well, as the meeting point
#' of two linearizations. At the minimum of \eqn{\lvert x \rvert} the bundle
#' holds subgradients near \eqn{+1} and near \eqn{-1}, and \eqn{0} is a
#' convex combination of them, which is the statement that \eqn{0} lies in
#' the subdifferential.
#'
#' The step solves
#' \deqn{\min_d\ \max_j \{ -\alpha_j + g_j'd \} + \frac{1}{2t}\lVert d \rVert^2,}
#' where \eqn{\alpha_j \ge 0} measures how far linearization \eqn{j} lies
#' below the objective at the current point, and the quadratic term keeps the
#' step inside the region where the model is trusted. The minimizer is
#' \eqn{d = -t\,p}, with \eqn{p} the aggregate subgradient and \eqn{\alpha}
#' the aggregate linearization error. The value of the piecewise-linear model
#' at \eqn{d}, \eqn{v = -(t\lVert p \rVert^2 + \alpha) \le 0}, is the
#' *predicted decrease*, which the acceptance test below reads.
#'
#' [crit_stationary()] reads a related but different quantity, the
#' *optimality estimate* \eqn{\lVert p \rVert^2 + \alpha}. Both vanish at a
#' solution, but \eqn{-v} contains the factor \eqn{t}, and \eqn{t} is halved
#' at every null step, so \eqn{-v} can fall below any tolerance because the
#' trust parameter has shrunk and not because the point is stationary. A run
#' that stopped on \eqn{-v} would then report convergence at a point where
#' the model still has a large aggregate subgradient. The optimality estimate
#' omits the factor \eqn{t} and does not depend on the size of the trust
#' parameter.
#'
#' # Serious steps and null steps
#'
#' A trial point is accepted, making a *serious step*, only if the decrease
#' that it delivers is at least the fraction `m_serious` of the predicted
#' decrease. This plays the role of the Armijo condition in a line search:
#' without a required fraction, a sequence of ever smaller improvements could
#' be accepted although it does not converge to a minimizer.
#'
#' A rejected trial point still contributes its subgradient to the bundle,
#' which improves the model at the next iteration. Such an iteration is a
#' *null step*, and a run with many null steps is refining the model of a
#' kink. Both counts appear in the `message` of the result, as in
#' `8 serious, 14 null`.
#'
#' `t` is halved after a null step and doubled after a serious one, within
#' `t_min` and `t_max`. Kiwiel's rule chooses the factor from a curvature
#' estimate and performs better; the rule implemented here is simpler and
#' bounded.
#'
#' # The trust parameter t0
#'
#' The step is \eqn{d = -t\,p}, so with a bare `t` the length of the first
#' step would be proportional to the size of the gradient. On a steep
#' objective that step lands where the objective is enormous, each null step
#' halves `t` while the gradients collected there are larger still, and the
#' run can spend its whole budget on rejected steps.
#'
#' When the largest absolute component of the gradient at the starting point
#' exceeds 1, `t0` is therefore divided by that component, so that the
#' largest coordinate of the first step has size `t0`; otherwise `t0` is used
#' as it is. This is the normalization of a line search that starts at 1
#' along a unit direction, and `t_min` and `t_max` are scaled in the same
#' way. A larger `t0` suits a problem whose optimum is far away, and a
#' smaller one a problem on which the model is accurate only nearby.
#'
#' Should a subgradient overflow anyway, which an objective that grows fast
#' enough can produce, the run stops with a message instead of spending its
#' budget on a subproblem whose matrix contains an infinity.
#'
#' # Bounded memory
#'
#' Without a limit the bundle grows indefinitely. When a new linearization
#' would make it exceed `bundle_size`, the oldest linearizations are replaced
#' by the *aggregate*, the affine function defined by the solution of the
#' subproblem at the current stability center. Discarding them would lose
#' their information and can stall the method; the aggregate retains a
#' summary of that information in one element.
#'
#' # Convexity requirement
#'
#' The theory is for convex \eqn{f}, where the linearization errors
#' \eqn{\alpha_j} are non-negative automatically. On a non-convex objective
#' they can be negative and are clipped at zero. The method still runs and
#' usually works, but the convergence guarantee no longer applies: a negative
#' error means that \eqn{f} lies below one of its own linearizations, and
#' setting the error to zero discards that information.
#'
#' # Subgradients
#'
#' Supply `gr`. Without it the package computes the gradient by finite
#' differences, which is accurate wherever \eqn{f} is differentiable, that
#' is, almost everywhere. A difference taken *across* a kink is not a
#' subgradient, and the model is then built from a number that does not
#' correspond to any linearization of \eqn{f}.
#'
#' @return An S7 object of class [Bundle], inheriting from [optimizer()], to
#'   be handed to [minimize()].
#'
#' @references
#' Kiwiel, K. C. (1990). Proximity control in bundle methods for convex
#' nondifferentiable minimization. *Mathematical Programming* **46**,
#' 105--122.
#'
#' Mäkelä, M. M. (2002). Survey of bundle methods for nonsmooth optimization.
#' *Optimization Methods and Software* **17**, 1--29.
#'
#' @examples
#' bundle()
#'
#' # The median is the minimizer of a sum of absolute deviations: the
#' # objective has a kink at every observation, and one of them is the answer.
#' set.seed(1)
#' y <- rnorm(101)
#' f <- function(p) sum(abs(y - p))
#' g <- function(p) -sum(sign(y - p))
#'
#' r <- minimize(bundle(), f, par = 0, gr = g)
#' c(bundle = r@par, median = median(y))
#' abs(r@par - median(y))     # 8e-16
#' r@message                  # how many trials were accepted
#'
#' # The aggregate subgradient goes to zero at the kink, whereas the
#' # subgradient computed by g just beside it has norm 1.
#' r@gradient
#' abs(g(r@par - 1e-9))
#'
#' # A descent method lands on the same point and cannot certify it.
#' d <- minimize(bfgs(criterion = crit_grad()), f, par = 0, gr = g)
#' c(par = d@par, converged = d@converged)
#'
#' # The kink lies along a diagonal. The bundle method and the simplex reach the
#' # minimum; the random-direction pattern search (the default of compass())
#' # comes close, whereas a coordinate search stops at a wrong point.
#' ff <- function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2)
#' gg <- function(p) c(sign(p[1] + p[2]), sign(p[1] + p[2])) + 0.2 * p
#' c(bundle      = minimize(bundle(), ff, c(1, 0.5), gr = gg)@value,
#'   nelder_mead = minimize(nelder_mead(), ff, c(1, 0.5))@value,
#'   compass     = minimize(compass(), ff, c(1, 0.5))@value)
#'
#' # A gradient rule is rejected, since this method declares no gradient.
#' try(minimize(bundle(criterion = crit_grad()), f, par = 0, gr = g))
#'
#' @seealso [prox_grad()] for the non-smooth case where a proximal operator is
#'   available, [nelder_mead()] and [compass()] for the derivative-free
#'   methods, [crit_stationary()] for the default rule of this method.
#' @export
bundle <- function(criterion = crit_stationary(),
                   t0 = 1, t_min = 1e-10, t_max = 1e10,
                   m_serious = 0.1, bundle_size = 20,
                   qp_iters = 500, qp_tol = 1e-12,
                   maxit = 500, max_eval = Inf,
                   verbose = FALSE, refresh = 10, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_tol(t0, "t0")
  check_tol(t_min, "t_min")
  check_tol(t_max, "t_max")
  check_tol(qp_tol, "qp_tol")
  if (t_min >= t_max) {
    stop("'t_min' must be strictly below 't_max'.", call. = FALSE)
  }
  if (t0 < t_min || t0 > t_max) {
    stop("'t0' must lie between 't_min' and 't_max'.", call. = FALSE)
  }
  if (length(m_serious) != 1L || !is.numeric(m_serious) || is.na(m_serious) ||
      m_serious <= 0 || m_serious >= 1) {
    stop("'m_serious' must be a single number in (0, 1).", call. = FALSE)
  }
  if (length(bundle_size) != 1L || !is.numeric(bundle_size) ||
      is.na(bundle_size) || bundle_size < 2 ||
      bundle_size > .Machine$integer.max) {
    stop("'bundle_size' must be a single finite number of at least 2.",
         call. = FALSE)
  }
  if (length(qp_iters) != 1L || !is.numeric(qp_iters) || is.na(qp_iters) ||
      qp_iters < 1) {
    stop("'qp_iters' must be a single positive number.", call. = FALSE)
  }
  Bundle(
    name = "proximal bundle", criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    t0 = t0, t_min = t_min, t_max = t_max, m_serious = m_serious,
    bundle_size = bundle_size, qp_iters = qp_iters, qp_tol = qp_tol
  )
}


#' @title What the Bundle Method Can Offer a Stopping Rule
#' @name optimizer_provides.Bundle
#'
#' @description
#' Reports `"stationarity"` and not `"gradient"`, although the method
#' evaluates subgradients and reports their aggregate. A rule built on
#' [crit_grad()] would test the subgradient at the iterate, which does not
#' tend to zero: near the minimum of \eqn{\lvert x \rvert} every evaluated
#' subgradient has norm 1, so the rule would not fire when the run reached
#' the minimizer.
#'
#' @details
#' [crit_stationary()] reads the optimality estimate
#' \eqn{\lVert p \rVert^{2} + \alpha} instead, which tends to zero, and is
#' the default rule of this method. [check_optimizer()] also reads this
#' declaration before its second check, which therefore does not require the
#' `gradient` reported by a bundle run to equal the gradient of the
#' objective.
#'
#' @param optimizer A `Bundle` object.
#'
#' @return The character vector `"stationarity"`.
#'
#' @examples
#' optimizer_provides(bundle())
#'
#' # A gradient rule is therefore rejected when the run starts.
#' f <- function(p) sum(abs(p - c(1, -2)))
#' try(minimize(bundle(criterion = crit_grad()), f, c(0, 0),
#'              gr = function(p) sign(p - c(1, -2))))
#'
#' @keywords internal
S7::method(optimizer_provides, Bundle) <- function(optimizer)
  "stationarity"


#' @title Minimize by the Proximal Bundle Method
#' @name minimize.Bundle
#'
#' @description
#' Runs [bundle()] on the objective: a cutting-plane model built from the
#' subgradients collected so far, a proximal subproblem solved for the step,
#' and an acceptance test that turns each trial into a serious step or a null
#' one.
#'
#' @param optimizer A `Bundle` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `gr` should return
#'   a subgradient, not a difference taken across a kink; `he` is
#'   accepted and ignored. Bounds are taken and removed by
#'   reparametrization.
#'
#' @return An [optimizer_result()] whose `gradient` is the aggregate
#'   subgradient, which tends to zero at a solution, and whose `message`
#'   reports the numbers of serious and null steps.
#'
#' @keywords internal
S7::method(minimize, Bundle) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    spec <- prepare_objective(optimizer, fn, par, gr, he)
    bounds <- check_bounds(lower, upper, par)

    t0 <- proc.time()[["elapsed"]]
    out <- bundle_run(
      spec = spec, par = as.numeric(par),
      criterion = optimizer@criterion, crit_fn = crit_met,
      t0 = optimizer@t0, t_min = optimizer@t_min, t_max = optimizer@t_max,
      m_serious = optimizer@m_serious,
      bundle_size = as.integer(optimizer@bundle_size),
      qp_iters = as.integer(optimizer@qp_iters), qp_tol = optimizer@qp_tol,
      maxit = as.integer(optimizer@maxit),
      max_eval = budget_int(optimizer@max_eval),
      verbose = optimizer@verbose,
      refresh = as.integer(optimizer@refresh),
      keep_trace = optimizer@keep_trace,
      bounds = bounds
    )
    build_result(out, optimizer, spec, proc.time()[["elapsed"]] - t0)
  }
