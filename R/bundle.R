#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

#' @title S7 Class for the Proximal Bundle Method
#' @description The class \code{\link{bundle}} instantiates.
#' @param t0 Initial proximity weight.
#' @param t_min,t_max Bounds on it.
#' @param m_serious Fraction of the predicted decrease a serious step must
#'   achieve.
#' @param bundle_size Largest number of linearisations kept.
#' @param qp_iters,qp_tol Effort spent on the subproblem.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{bundle}}
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
#' The method for a convex objective with kinks in it: it uses subgradients,
#' and unlike every descent method here it is not misled by them.
#'
#' @param criterion The stopping rule. Defaults to
#'   \code{crit_stationary(1e-8)}, on the predicted decrease.
#' @param t0 Initial proximity weight, expressed as a \strong{step length} on
#'   the parameter scale rather than as a bare multiplier; see Details.
#'   Defaults to \code{1}.
#' @param t_min,t_max Bounds on it, so that neither a run of null steps nor a
#'   run of serious ones can drive it to zero or to infinity. Default
#'   \code{1e-10} and \code{1e10}.
#' @param m_serious Fraction of the predicted decrease that a step must actually
#'   deliver to be accepted, in \eqn{(0, 1)}. Defaults to \code{0.1}.
#' @param bundle_size Largest number of linearisations kept before the oldest
#'   are replaced by their aggregate. Defaults to \code{20}.
#' @param qp_iters,qp_tol Effort spent on the subproblem: at most this many
#'   accelerated projected-gradient steps, stopping when the weights move by
#'   less than \code{qp_tol}. Defaults \code{500} and \code{1e-12}.
#' @param maxit Maximum iterations. Defaults to 500.
#' @param max_eval Maximum objective evaluations. Defaults to 10000.
#' @param verbose Report progress? Defaults to \code{FALSE}.
#' @param refresh Report every this many iterations. Defaults to 10.
#' @param keep_trace Store the iteration path? Defaults to \code{FALSE}.
#'
#' @details
#' \subsection{Behaviour at a kink}{
#' At the minimum of \eqn{\lvert x \rvert} the evaluated subgradient is
#' \eqn{\pm 1}. A descent method therefore sees a large gradient, proposes a
#' step, finds no acceptable one, and stops reporting failure while standing
#' exactly on the answer — and no tolerance can fix that, because the quantity
#' it is testing does not become small.
#'
#' The bundle method does not test any single subgradient. It keeps a
#' \emph{collection} of them, from the points it has visited, and builds the
#' piecewise-linear model \eqn{\max_j \{ f_j + g_j'(x - x_j) \}} — the largest of
#' the linearisations. A kink is exactly what such a model represents well: two
#' linearisations meeting. At the minimum of \eqn{\lvert x \rvert} the bundle
#' holds subgradients near \eqn{+1} and near \eqn{-1}, and \eqn{0} is a convex
#' combination of them, which is precisely the statement that \eqn{0} lies in
#' the subdifferential.
#'
#' That is what is tested. The step solves
#' \deqn{\min_d\ \max_j \{ -\alpha_j + g_j'd \} + \frac{1}{2t}\lVert d \rVert^2,}
#' where \eqn{\alpha_j \ge 0} measures how badly linearisation \eqn{j} misses
#' the current point and the quadratic term keeps the step inside the region
#' where the model is believed. Its optimal value \eqn{v \le 0} is the
#' \emph{predicted decrease}, which is what the acceptance test below uses.
#'
#' What \code{\link{crit_stationary}} watches is a related but different
#' quantity, the \emph{optimality estimate}
#' \eqn{\lVert p \rVert^2 + \alpha}, where \eqn{p} is the aggregate
#' subgradient. Both vanish together at a solution, and the distinction is not
#' pedantry: \eqn{-v} carries a factor of \eqn{t}, and \eqn{t} is halved at
#' every null step, so it can be pushed below any tolerance by the trust
#' parameter shrinking rather than by the point becoming stationary. A run doing
#' that reports success while standing somewhere its own model still calls
#' steep. Dropping the factor leaves the model's own claim, which no amount of
#' shrinking can flatter.
#' }
#'
#' \subsection{Serious steps and null steps}{
#' A trial point is accepted — a \emph{serious step} — only if it delivers at
#' least \code{m_serious} of the decrease the model promised. That is the same
#' bargain an Armijo condition strikes, for the same reason: accepting any
#' decrease at all lets a sequence of ever tinier improvements masquerade as
#' progress.
#'
#' When it is refused, the iteration is not wasted. The trial point contributes
#' its subgradient to the bundle, so the model is strictly better next time;
#' this is a \emph{null step}, and a run that reports many of them is refining
#' its picture of a kink, not failing. Both counts appear in the result's
#' message.
#'
#' \code{t} is halved after a null step and doubled after a serious one, within
#' \code{t_min} and \code{t_max}. Kiwiel's rule chooses the factor from a
#' curvature estimate and is better; this is the safeguarded version, which is
#' cruder and bounded.
#' }
#'
#' \subsection{The trust parameter t0}{
#' The step is \eqn{d = -t\,p}, so a bare \code{t} would make the first step as
#' long as the gradient happens to be big. On Rosenbrock from its usual start
#' that is 233 units, landing where the objective is \eqn{10^{11}} and its
#' gradient \eqn{10^{15}}; the null step that follows halves \code{t} while the
#' gradient has squared, so \code{t} can never catch up and the run spends its
#' entire budget on rejected steps before returning the point it began at.
#'
#' \code{t0} is therefore divided by the size of the gradient at the starting
#' point, which makes the first step of length \code{t0} in parameter space.
#' It is the same normalisation a line search performs when it starts at 1 along
#' a unit direction, and \code{t_min} and \code{t_max} move with it since they
#' bound the same quantity. Raise \code{t0} for a problem whose optimum is far
#' away, lower it for one where the model is trustworthy only nearby.
#'
#' Should a subgradient overflow anyway — possible for an objective that grows
#' fast enough — the run stops and says so, rather than spending its budget on a
#' subproblem whose matrix contains an infinity.
#' }
#'
#' \subsection{Bounded memory}{
#' Left alone the bundle grows without limit. When it reaches
#' \code{bundle_size} the oldest linearisations are replaced by the
#' \emph{aggregate} — the single affine function the subproblem's solution
#' defines — rather than discarded. Discarding loses what they knew and can stall
#' the method; the aggregate keeps a summary of all of it in one element, which
#' is what makes a bounded bundle safe.
#' }
#'
#' \subsection{Convexity requirement}{
#' The theory is for convex \eqn{f}, where the linearisation errors
#' \eqn{\alpha_j} are non-negative automatically. On a non-convex objective they
#' can come out negative and are clipped at zero. The method then still runs and
#' usually works, but that clip is exactly where the guarantee stops: a negative
#' error is the model reporting that \eqn{f} lies below its own linearisation,
#' and setting it to zero suppresses the information rather than using it.
#' }
#'
#' \subsection{Subgradients}{
#' Supply \code{gr}. Without one the package differences the objective, and
#' although a difference is a perfectly good gradient wherever \eqn{f} is
#' differentiable — which is almost everywhere, so a run mostly gets away with
#' it — a difference taken \emph{across} a kink is not a subgradient of anything
#' and the model will be built from a number that belongs to no linearisation.
#' }
#'
#' @return An S7 object of class \code{Bundle}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @references
#' Kiwiel, K. C. (1990). Proximity control in bundle methods for convex
#' nondifferentiable minimization. \emph{Mathematical Programming} \strong{46},
#' 105--122.
#'
#' Mäkelä, M. M. (2002). Survey of bundle methods for nonsmooth optimization.
#' \emph{Optimization Methods and Software} \strong{17}, 1--29.
#'
#' @examples
#' bundle()
#'
#' # the median, as the minimiser of a sum of absolute deviations: the objective
#' # has a kink at every observation, and one of them is the answer
#' set.seed(1)
#' y <- rnorm(101)
#' f <- function(p) sum(abs(y - p))
#' g <- function(p) -sum(sign(y - p))
#' r <- minimize(bundle(), f, par = 0, gr = g)
#' c(bundle = r@par, median = median(y))
#'
#' # a kink running diagonally, where a coordinate-wise search stalls
#' minimize(bundle(), function(p) abs(p[1] + p[2]) + 0.1 * sum(p^2),
#'          c(1, 0.5), gr = function(p) c(sign(p[1] + p[2]), sign(p[1] + p[2])) +
#'                          0.2 * p)@value
#'
#' @seealso \code{\link{nelder_mead}}, \code{\link{compass}},
#'   \code{\link{crit_stationary}}
#' @export
bundle <- function(criterion = crit_stationary(1e-8),
                   t0 = 1, t_min = 1e-10, t_max = 1e10,
                   m_serious = 0.1, bundle_size = 20,
                   qp_iters = 500, qp_tol = 1e-12,
                   maxit = 500, max_eval = 10000,
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
#' @description
#' The objective and the predicted decrease, but not a gradient.
#' @details
#' It evaluates subgradients and reports their aggregate, but it does not offer
#' \code{"gradient"}, and the omission is deliberate. \code{\link{crit_grad}}
#' would then test a quantity that does not go to zero — at the minimum of
#' \eqn{\lvert x \rvert} every subgradient has norm 1 — so the rule would never
#' fire while sitting on the answer. \code{\link{crit_stationary}} tests the
#' predicted decrease instead, which does go to zero.
#' @param optimizer A \code{Bundle} object.
#' @return A character vector.
#' @keywords internal
S7::method(optimizer_provides, Bundle) <- function(optimizer)
  "stationarity"


#' @title Minimise by the Proximal Bundle Method
#' @name minimize.Bundle
#' @description Runs \code{\link{bundle}} on the objective.
#' @param optimizer A \code{Bundle} object.
#' @param fn,par,gr,he,lower,upper,... As in \code{\link{minimize}}. \code{gr} should
#'   return a subgradient; \code{he} is ignored.
#' @return An \code{\link{optimizer_result}}, whose \code{gradient} is the
#'   \emph{aggregate} subgradient, the one quantity of that shape which goes to
#'   zero at a solution.
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
      max_eval = as.integer(optimizer@max_eval),
      verbose = optimizer@verbose,
      refresh = as.integer(optimizer@refresh),
      keep_trace = optimizer@keep_trace,
      bounds = bounds
    )
    build_result(out, optimizer, spec, proc.time()[["elapsed"]] - t0)
  }
