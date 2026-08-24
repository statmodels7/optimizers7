#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

#' @title S7 Class for the Proximal Bundle Method
#'
#' @description
#' An optimizer holding the proximity weight and its bounds, the acceptance
#' fraction that separates a serious step from a null one, the size of the
#' bundle of linearizations, and the effort spent on the subproblem. Built by
#' [bundle()].
#'
#' @details
#' Beyond the seven properties every optimizer has, a `Bundle` carries seven
#' of its own, in three groups: the proximity weight (`t0`, `t_min`,
#' `t_max`), the acceptance test (`m_serious`), and the model and its
#' subproblem (`bundle_size`, `qp_iters`, `qp_tol`).
#'
#' @param t0 Initial proximity weight, as a step length on the parameter
#'   scale.
#' @param t_min,t_max Bounds on the proximity weight.
#' @param m_serious Fraction of the predicted decrease a serious step must
#'   achieve.
#' @param bundle_size Largest number of linearizations kept.
#' @param qp_iters,qp_tol Effort spent on the subproblem.
#'
#' @return An S7 object of class `Bundle` inheriting from [optimizer()], with
#'   the seven properties above beside the seven shared ones.
#'
#' @seealso [bundle()] for the constructor, [crit_stationary()] for the rule
#'   it reads.
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
#' term keeps the step near the stability center, and convergence is tested on
#' the aggregate subgradient rather than on any single one, which does not
#' vanish at a kink.
#'
#' @param criterion The stopping rule. Defaults to
#'   `crit_stationary(1e-8)`, on the predicted decrease.
#' @param t0 Initial proximity weight, expressed as a **step length** on
#'   the parameter scale rather than as a bare multiplier; see Details.
#'   Defaults to `1`.
#' @param t_min,t_max Bounds on it, so that neither a run of null steps nor a
#'   run of serious ones can drive it to zero or to infinity. Default
#'   `1e-10` and `1e10`.
#' @param m_serious Fraction of the predicted decrease that a step must actually
#'   deliver to be accepted, in \eqn{(0, 1)}. Defaults to `0.1`.
#' @param bundle_size Largest number of linearizations kept before the oldest
#'   are replaced by their aggregate. Defaults to `20`.
#' @param qp_iters,qp_tol Effort spent on the subproblem: at most this many
#'   accelerated projected-gradient steps, stopping when the weights move by
#'   less than `qp_tol`. Defaults `500` and `1e-12`.
#' @param maxit Maximum iterations. Defaults to 500.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`:
#'   no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 10.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#'
#' @details
#' # Behavior at a kink
#'
#' At the minimum of \eqn{\lvert x \rvert} the evaluated subgradient is
#' \eqn{\pm 1}. A descent method therefore sees a large gradient, proposes a
#' step, finds no acceptable one, and stops while standing exactly on the
#' answer. No tolerance can fix that: the quantity being tested never becomes
#' small. Measured on the median example below, [bfgs()] and [lbfgs()] both
#' land on the median to six figures and report `converged = FALSE` with
#' `stopped without converging`.
#'
#' The bundle method tests no single subgradient. It keeps a *collection* of
#' them, from the points it has visited, and builds the piecewise-linear
#' model \eqn{\max_j \{ f_j + g_j'(x - x_j) \}}, the largest of the
#' linearizations. A kink is exactly what such a model represents well: two
#' linearizations meeting. At the minimum of \eqn{\lvert x \rvert} the bundle
#' holds subgradients near \eqn{+1} and near \eqn{-1}, and \eqn{0} is a
#' convex combination of them, which is precisely the statement that \eqn{0}
#' lies in the subdifferential.
#'
#' That is what is tested. The step solves
#' \deqn{\min_d\ \max_j \{ -\alpha_j + g_j'd \} + \frac{1}{2t}\lVert d \rVert^2,}
#' where \eqn{\alpha_j \ge 0} measures how badly linearization \eqn{j} misses
#' the current point and the quadratic term keeps the step inside the region
#' where the model is believed. Its optimal value \eqn{v \le 0} is the
#' *predicted decrease*, which is what the acceptance test below uses.
#'
#' What [crit_stationary()] watches is a related but different
#' quantity, the *optimality estimate*
#' \eqn{\lVert p \rVert^2 + \alpha}, where \eqn{p} is the aggregate
#' subgradient. Both vanish together at a solution, and the distinction is not
#' pedantry: \eqn{-v} carries a factor of \eqn{t}, and \eqn{t} is halved at
#' every null step, so it can be pushed below any tolerance by the trust
#' parameter shrinking rather than by the point becoming stationary. A run doing
#' that reports success while standing somewhere its own model still calls
#' steep. Dropping the factor leaves the model's own claim, which no amount of
#' shrinking can flatter.
#'
#' # Serious steps and null steps
#'
#' A trial point is accepted, making a *serious step*, only if it delivers at
#' least `m_serious` of the decrease the model promised. That is the same
#' bargain an Armijo condition strikes, for the same reason: accepting any
#' decrease at all lets a sequence of ever tinier improvements masquerade as
#' progress.
#'
#' A rejected trial is not a wasted iteration. The trial point contributes its
#' subgradient to the bundle, so the model is strictly better next time. That
#' is a *null step*, and a run reporting many of them is refining its picture
#' of a kink. Both counts appear in the result's `message`, as in
#' `8 serious, 14 null`.
#'
#' `t` is halved after a null step and doubled after a serious one, within
#' `t_min` and `t_max`. Kiwiel's rule chooses the factor from a curvature
#' estimate and is better; this is the safeguarded version, cruder and
#' bounded.
#'
#' # The trust parameter t0
#' The step is \eqn{d = -t\,p}, so a bare `t` would make the first step as
#' long as the gradient happens to be big. On Rosenbrock from its usual start
#' the gradient has norm 232.9, and a step of that length lands where the
#' objective is \eqn{2.1 \times 10^{11}}. The null step that follows halves
#' `t` while the gradient has squared, so `t` can never catch up and the run
#' spends its whole budget on rejected steps before returning the point it
#' began at.
#'
#' `t0` is therefore divided by the size of the gradient at the starting
#' point, which makes the first step of length `t0` in parameter space. It is
#' the same normalization a line search performs when it starts at 1 along a
#' unit direction, and `t_min` and `t_max` move with it, bounding the same
#' quantity. Raise `t0` for a problem whose optimum is far away, lower it for
#' one where the model is trustworthy only nearby.
#'
#' Should a subgradient overflow anyway, which an objective that grows fast
#' enough can produce, the run stops with a message instead of spending its
#' budget on a subproblem whose matrix contains an infinity.
#'
#' # Bounded memory
#'
#' Left alone the bundle grows without limit. When it reaches `bundle_size`
#' the oldest linearizations are replaced by the *aggregate*, the single
#' affine function the subproblem's solution defines. Discarding them instead
#' would lose what they knew and can stall the method; the aggregate keeps a
#' summary of all of it in one element, and that is what makes a bounded
#' bundle safe.
#'
#' # Convexity requirement
#'
#' The theory is for convex \eqn{f}, where the linearization errors
#' \eqn{\alpha_j} are non-negative automatically. On a non-convex objective
#' they can come out negative and are clipped at zero. The method still runs
#' and usually works, and that clip is exactly where the guarantee stops: a
#' negative error is the model reporting that \eqn{f} lies below its own
#' linearization, and setting it to zero suppresses the information instead of
#' using it.
#'
#' # Subgradients
#'
#' Supply `gr`. Without one the package differences the objective, and a
#' difference is a perfectly good gradient wherever \eqn{f} is
#' differentiable, which is almost everywhere, so a run mostly gets away with
#' it. A difference taken *across* a kink is a subgradient of nothing, and the
#' model is then built from a number that belongs to no linearization.
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
#' # The median, as the minimizer of a sum of absolute deviations: the
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
#' # The aggregate subgradient does go to zero at the kink, where every single
#' # subgradient the objective offers there has norm 1.
#' r@gradient
#' abs(g(r@par - 1e-9))
#'
#' # A descent method lands on the same point and cannot certify it.
#' d <- minimize(bfgs(criterion = crit_grad()), f, par = 0, gr = g)
#' c(par = d@par, converged = d@converged)
#'
#' # A kink running diagonally. The bundle method and the simplex both find it;
#' # a coordinatewise search cannot, its axes never crossing the kink.
#' ff <- function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2)
#' gg <- function(p) c(sign(p[1] + p[2]), sign(p[1] + p[2])) + 0.2 * p
#' c(bundle      = minimize(bundle(), ff, c(1, 0.5), gr = gg)@value,
#'   nelder_mead = minimize(nelder_mead(), ff, c(1, 0.5))@value,
#'   compass     = minimize(compass(), ff, c(1, 0.5))@value)
#'
#' # A gradient rule is refused, this method offering no gradient.
#' try(minimize(bundle(criterion = crit_grad()), f, par = 0, gr = g))
#'
#' @seealso [prox_grad()] for the non-smooth case where a proximal operator is
#'   available, [nelder_mead()] and [compass()] for the derivative-free
#'   methods, [crit_stationary()] for the rule this one reads.
#' @export
bundle <- function(criterion = crit_stationary(),
                   t0 = 1, t_min = 1e-10, t_max = 1e10,
                   m_serious = 0.1, bundle_size = 20,
                   qp_iters = 500, qp_tol = 1e-12,
                   maxit = 500, max_eval = Inf,
                   verbose = FALSE, refresh = 10, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)
  check_tol(t0)
  check_tol(t_min)
  check_tol(t_max)
  check_tol(qp_tol)
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
      is.na(bundle_size) || bundle_size < 2) {
    stop("'bundle_size' must be a single number of at least 2.", call. = FALSE)
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
#' Reports `"stationarity"` and withholds `"gradient"`, though the method
#' evaluates subgradients and reports their aggregate. The omission is
#' deliberate: [crit_grad()] would test a quantity that never goes to zero,
#' since at the minimum of \eqn{\lvert x \rvert} every subgradient has norm 1,
#' so the rule would sit on the answer without firing.
#'
#' @details
#' [crit_stationary()] reads the optimality estimate
#' \eqn{\lVert p \rVert^{2} + \alpha} instead, which does go to zero, and is
#' this method's default rule. [check_optimizer()] also consults this
#' function before its second check, so a bundle run is not asked to prove
#' that what it reports as `gradient` is one.
#'
#' @param optimizer A `Bundle` object.
#'
#' @return The character vector `"stationarity"`.
#'
#' @examples
#' optimizer_provides(bundle())
#'
#' # So a gradient rule is refused when the run starts.
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
#'   a **subgradient**, not a difference taken across a kink; `he` is
#'   accepted and ignored. Bounds are taken and removed by
#'   reparametrization.
#'
#' @return An [optimizer_result()] whose `gradient` is the **aggregate**
#'   subgradient, the one quantity of that shape which goes to zero at a
#'   solution, and whose `message` reports the serious and null step counts.
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
