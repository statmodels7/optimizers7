#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

#' @title S7 Class for Adam
#' @description The class \code{\link{adam}} instantiates.
#' @param alpha The learning rate.
#' @param beta1,beta2 Decay rates for the two moment estimates.
#' @param eps The denominator floor.
#' @param decay Rate at which the learning rate is reduced.
#' @param amsgrad Whether to hold the second moment at its running maximum.
#' @param resample Fraction of the observations used at each iteration.
#' @return An S7 object inheriting from \code{\link{optimizer}}.
#' @seealso \code{\link{adam}}
#' @keywords internal
Adam <- S7::new_class("Adam", parent = optimizer,
  properties = list(
    alpha    = S7::class_numeric,
    beta1    = S7::class_numeric,
    beta2    = S7::class_numeric,
    eps      = S7::class_numeric,
    decay    = S7::class_numeric,
    amsgrad  = S7::class_logical,
    resample = S7::class_numeric
  ))


#' @title Adaptive Moment Estimation
#'
#' @description
#' Adam: a first-order method that gives every coordinate its own step length,
#' inferred from the size of the gradients it has been seeing. Available both on
#' the whole sample and on subsamples of it.
#'
#' @param criterion The stopping rule. Defaults to \code{\link{crit_never}}, so
#'   the run is governed by \code{maxit}; see Details.
#' @param alpha The learning rate — the size of a step when the gradient is
#'   steady. Defaults to \code{0.01}.
#' @param beta1 Decay rate of the first moment, the smoothed gradient. Defaults
#'   to \code{0.9}.
#' @param beta2 Decay rate of the second moment, the smoothed squared gradient.
#'   Defaults to \code{0.999}.
#' @param eps Added to the square-rooted second moment before dividing, so that
#'   a coordinate whose gradient has been uniformly zero does not divide by it.
#'   Defaults to \code{1e-8}.
#' @param decay Reduces the learning rate as \eqn{\alpha_t = \alpha/(1 + \delta
#'   t)}. Defaults to \code{0}, a constant rate; see Details.
#' @param amsgrad Hold the second moment at its running maximum? Defaults to
#'   \code{FALSE}; see Details.
#' @param resample Fraction of the terms to draw at each iteration, in
#'   \code{(0, 1]}. Defaults to \code{1}, the whole sample. Anything less
#'   requires a \code{\link{finite_sum}} objective.
#' @param maxit Maximum iterations. Defaults to 1000, higher than the other
#'   methods because Adam takes many small steps rather than few large ones.
#' @param max_eval Maximum objective evaluations. Defaults to 100000.
#' @param verbose Report progress? Defaults to \code{FALSE}.
#' @param refresh Report every this many iterations. Defaults to 100.
#' @param keep_trace Store the iteration path? Defaults to \code{FALSE}.
#'
#' @details
#' The idea is one line. Adam keeps an exponentially weighted average of the
#' gradient, \eqn{m_t}, and of its square, \eqn{v_t}, and steps
#' \deqn{x_{t+1} = x_t - \alpha\, \hat m_t / (\sqrt{\hat v_t} + \epsilon).}
#' The division is elementwise, which is the whole point: a coordinate whose
#' gradient has been consistently large is divided by a large number and moves
#' modestly, while one whose gradient is small but persistent still moves. It is
#' a diagonal preconditioner assembled from the gradients already seen, so it
#' costs nothing beyond them — and that is also its limit, since a diagonal
#' cannot represent the correlation between parameters that \code{\link{bfgs}}
#' learns from the same information.
#'
#' Both averages start at zero, so early on they are pulled towards it; dividing
#' by \eqn{1 - \beta^t} removes exactly that bias, and without it the first
#' iterations would barely move.
#'
#' \subsection{Why it is not a descent method}{
#' Adam takes no line search and makes no attempt to decrease the objective at
#' every step. That is deliberate, not an omission: the freedom to go uphill is
#' most of why it tolerates a gradient computed from a handful of observations.
#' It also means that none of the usual reassurances apply. There is no
#' guarantee of monotone progress, and the run may end somewhere worse than it
#' passed through.
#'
#' The practical consequence is that Adam is the wrong tool for a small smooth
#' problem where a Hessian is affordable. Use \code{\link{newton}} or
#' \code{\link{bfgs}} there and reach machine precision in a dozen iterations.
#' Adam earns its place when the parameter vector is long, when the objective is
#' a sum over many observations, or when the surface is rough enough that a
#' quadratic model is a fiction.
#' }
#'
#' \subsection{Stopping, and why the default is a budget}{
#' The default criterion is \code{\link{crit_never}}: the run ends when
#' \code{maxit} is reached, and reports \code{converged = FALSE}, which is the
#' truth about a run that nothing checked.
#'
#' With \code{resample = 1} the objective and the gradient are exact, and a real
#' rule such as \code{crit_grad(1e-6)} may be passed and will work. With
#' \code{resample < 1} both are estimates that depend on which observations were
#' drawn, so a tolerance applied to either would be measuring the sampling noise
#' rather than the progress. Such a criterion is therefore \strong{refused},
#' by name, rather than accepted and left never to fire — the same discipline as
#' everywhere else in the toolkit. Rules on the movement of the parameters
#' remain available, and with \code{decay > 0} they are meaningful.
#' }
#'
#' \subsection{The safeguards}{
#' \code{eps} floors the denominator. \code{decay} makes the learning rate
#' \eqn{O(1/t)}, which is the Robbins–Monro condition a stochastic run needs to
#' settle at the optimum rather than rattle about it at a radius set by
#' \eqn{\alpha}; it is off by default because with \code{resample = 1} there is
#' nothing to average away.
#'
#' \code{amsgrad} replaces \eqn{v_t} by its running maximum. This is not
#' cosmetic: Reddi, Kale and Kumar (2018) exhibited a convex problem on which
#' Adam as published fails to converge, because \eqn{v_t} can shrink and let a
#' single large gradient dominate the iterate long after it has passed. The
#' maximum forbids that, at the cost of steps that only ever get shorter. It is
#' \code{FALSE} by default so that \code{adam()} is Adam — a run that silently
#' did something else would not reproduce anything — but it is worth turning on
#' whenever a run refuses to settle.
#'
#' A non-finite gradient or update stops the run and says so, rather than
#' propagating a \code{NaN} into every iterate after it.
#' }
#'
#' \subsection{Subsampling}{
#' \code{resample} is the fraction of terms drawn at each iteration, and it
#' requires a \code{\link{finite_sum}} objective, which is what tells the
#' optimiser that the objective \emph{has} terms. A plain \code{fn(par)} is a
#' black box: the optimiser does not know what an observation is and must not
#' guess. An analytic gradient is required too, since differencing a subsample
#' would cost the \eqn{2p} evaluations that subsampling exists to avoid.
#'
#' Draws go through \R's generator, so \code{set.seed()} governs the run and it
#' can be reproduced. The value and gradient finally reported are recomputed on
#' the \strong{whole} sample at the point actually reached: two minibatch values
#' differ by which terms were drawn as much as by where the parameters went, so
#' anything else would let a lucky draw decide what the run claims.
#' }
#'
#' @return An S7 object of class \code{Adam}, inheriting from
#'   \code{\link{optimizer}}.
#'
#' @references
#' Kingma, D. P. and Ba, J. (2015). Adam: A Method for Stochastic Optimization.
#' \emph{ICLR}.
#'
#' Reddi, S. J., Kale, S. and Kumar, S. (2018). On the Convergence of Adam and
#' Beyond. \emph{ICLR}.
#'
#' @examples
#' adam()
#' adam(alpha = 0.05, amsgrad = TRUE)
#'
#' # full sample, on a quadratic
#' minimize(adam(alpha = 0.1, maxit = 2000),
#'          function(p) sum((p - c(1, 2))^2), c(0, 0),
#'          gr = function(p) 2 * (p - c(1, 2)))
#'
#' # minibatches, which need an objective that knows it is a sum
#' set.seed(1)
#' y <- rnorm(500, mean = 2)
#' obj <- finite_sum(fn = function(par, idx) sum((y[idx] - par)^2) / 2,
#'                   gr = function(par, idx) -sum(y[idx] - par),
#'                   n  = length(y))
#' minimize(adam(alpha = 0.05, resample = 0.1, decay = 0.01), obj, par = 0)
#'
#' @seealso \code{\link{bfgs}}, \code{\link{finite_sum}},
#'   \code{\link{crit_never}}
#' @export
adam <- function(criterion = crit_never(),
                 alpha = 0.01, beta1 = 0.9, beta2 = 0.999, eps = 1e-8,
                 decay = 0, amsgrad = FALSE, resample = 1,
                 maxit = 1000, max_eval = 100000,
                 verbose = FALSE, refresh = 100, keep_trace = FALSE) {
  check_optimizer_args(criterion, maxit, max_eval, verbose, refresh, keep_trace)

  pos <- function(v, nm) {
    if (length(v) != 1L || !is.numeric(v) || is.na(v) || v <= 0) {
      stop("'", nm, "' must be a single positive number.", call. = FALSE)
    }
  }
  pos(alpha, "alpha")
  pos(eps, "eps")
  unit <- function(v, nm) {
    if (length(v) != 1L || !is.numeric(v) || is.na(v) || v < 0 || v >= 1) {
      stop("'", nm, "' must be a single number in [0, 1).", call. = FALSE)
    }
  }
  unit(beta1, "beta1")
  unit(beta2, "beta2")
  if (length(decay) != 1L || !is.numeric(decay) || is.na(decay) || decay < 0) {
    stop("'decay' must be a single non-negative number.", call. = FALSE)
  }
  if (length(amsgrad) != 1L || !is.logical(amsgrad) || is.na(amsgrad)) {
    stop("'amsgrad' must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(resample) != 1L || !is.numeric(resample) || is.na(resample) ||
      resample <= 0 || resample > 1) {
    stop("'resample' must be a single number in (0, 1].", call. = FALSE)
  }

  Adam(
    name = "adam", criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    alpha = alpha, beta1 = beta1, beta2 = beta2, eps = eps,
    decay = decay, amsgrad = amsgrad, resample = resample
  )
}


#' @title What Adam Can Offer a Stopping Rule
#' @name optimizer_provides.Adam
#' @description
#' Everything, on the whole sample; nothing, on subsamples.
#' @details
#' With \code{resample < 1} both the objective and the gradient are estimates
#' from whichever terms were drawn. A tolerance on either measures the sampling
#' noise, so every rule that reads one is refused rather than accepted and left
#' to never fire.
#' @param optimizer An \code{Adam} object.
#' @return A character vector.
#' @keywords internal
S7::method(optimizer_provides, Adam) <- function(optimizer) {
  if (optimizer@resample < 1) character() else c("gradient", "objective")
}


#' @title Minimise by Adam
#' @name minimize.Adam
#' @description
#' Runs \code{\link{adam}} on the objective.
#' @param optimizer An \code{Adam} object.
#' @param fn,par,gr,he,bounds,... As in \code{\link{minimize}}.
#' @return An \code{\link{optimizer_result}}.
#' @keywords internal
S7::method(minimize, Adam) <-
  function(optimizer, fn, par, gr = NULL, he = NULL, bounds = NULL, ...) {
    stochastic <- optimizer@resample < 1

    spec <- prepare_objective(
      optimizer, fn, par, gr, he,
      note = if (stochastic) paste0(
        "With resample < 1 both are minibatch estimates, so a tolerance on ",
        "either would measure the sampling noise. Use crit_never() and a ",
        "budget, or a rule on the parameters.") else NULL
    )

    if (stochastic) {
      if (!identical(spec$kind, "finite_sum")) {
        stop("A resample below 1 needs a finite_sum() objective: a plain ",
             "function has no terms to draw from.", call. = FALSE)
      }
      if (!isTRUE(spec$has_gradient)) {
        stop("A resample below 1 needs an analytic gradient in finite_sum(): ",
             "differencing a subsample would cost the evaluations that ",
             "subsampling exists to save.", call. = FALSE)
      }
    }

    bounds <- check_bounds(bounds, par)

    # A subsampling run draws its minibatches from R's generator, so the state
    # it began with is the only thing that makes it repeatable.
    seed <- if (stochastic) capture_seed() else NULL

    # The value is not needed by the algorithm at all -- Adam steps on the
    # gradient alone -- so it is computed only when something will read it.
    need_value <- optimizer@keep_trace || optimizer@verbose ||
      "objective" %in% crit_needs(optimizer@criterion)

    t0 <- proc.time()[["elapsed"]]
    out <- adam_run(
      spec = spec, par = as.numeric(par),
      criterion = optimizer@criterion, crit_fn = crit_met,
      alpha = optimizer@alpha, beta1 = optimizer@beta1,
      beta2 = optimizer@beta2, eps = optimizer@eps,
      decay = optimizer@decay, amsgrad = optimizer@amsgrad,
      resample = optimizer@resample,
      maxit = as.integer(optimizer@maxit),
      max_eval = as.integer(optimizer@max_eval),
      verbose = optimizer@verbose,
      refresh = as.integer(optimizer@refresh),
      keep_trace = optimizer@keep_trace,
      need_value = need_value,
      bounds = bounds
    )
    elapsed <- proc.time()[["elapsed"]] - t0

    build_result(out, optimizer, spec, elapsed, seed)
  }
