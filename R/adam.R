#' @include optimizer_class.R
#' @include generics.R
#' @include methods.R
NULL

#' @title S7 Class for Adam
#'
#' @description
#' An optimizer holding the learning rate, the two moment decay rates and the
#' three safeguards of the Adam iteration. Built by [adam()]. It is the one
#' shipped method whose default stopping rule is [crit_never()], so a run
#' ends on its iteration budget and reports `converged = FALSE`.
#'
#' @details
#' Beyond the seven properties every optimizer has, an `Adam` carries six of
#' its own: `alpha`, `beta1` and `beta2` for the iteration as Kingma and Ba
#' published it, and `eps`, `decay` and `amsgrad` for the three repairs.
#'
#' @param alpha The learning rate, the size of a step when the gradient is
#'   steady.
#' @param beta1,beta2 Decay rates for the first and second moment estimates.
#' @param eps Added to the square-rooted second moment before dividing.
#' @param decay Rate at which the learning rate is reduced.
#' @param amsgrad Logical; whether the second moment is held at its running
#'   maximum.
#'
#' @return An S7 object of class `Adam` inheriting from [optimizer()], with
#'   the six properties above beside the seven shared ones.
#'
#' @seealso [adam()] for the constructor, [crit_never()] for its default
#'   rule.
#' @name Adam-class
#' @aliases Adam
#' @keywords internal
Adam <- S7::new_class("Adam", parent = optimizer,
  properties = list(
    alpha    = S7::class_numeric,
    beta1    = S7::class_numeric,
    beta2    = S7::class_numeric,
    eps      = S7::class_numeric,
    decay    = S7::class_numeric,
    amsgrad  = S7::class_logical
  ))


#' @title Adaptive Moment Estimation
#'
#' @description
#' Adaptive moment estimation: a first-order method whose coordinate-wise step
#' lengths come from exponentially weighted first and second moments of the
#' gradient. It takes no line search and tolerates a gradient that is correct
#' only on average, which makes it the method suited to a stochastic
#' objective.
#'
#' @param criterion The stopping rule. Defaults to [crit_never()], so
#'   the run is governed by `maxit`; see Details.
#' @param alpha The learning rate, a single positive number: the size of a
#'   step when the gradient is steady. Defaults to `0.01`.
#' @param beta1 Decay rate of the first moment, the smoothed gradient. Defaults
#'   to `0.9`.
#' @param beta2 Decay rate of the second moment, the smoothed squared gradient.
#'   Defaults to `0.999`.
#' @param eps Added to the square-rooted second moment before dividing, so that
#'   a coordinate whose gradient has been uniformly zero does not divide by it.
#'   Defaults to `1e-8`.
#' @param decay Reduces the learning rate as \eqn{\alpha_t = \alpha/(1 + \delta
#'   t)}. Defaults to `0`, a constant rate; see Details.
#' @param amsgrad Hold the second moment at its running maximum? Defaults to
#'   `FALSE`; see Details.
#' @param maxit Maximum iterations. Defaults to 1000, higher than the other
#'   methods because Adam takes many small steps where a second-order method
#'   takes few large ones.
#' @param max_eval Maximum objective evaluations. Defaults to `Inf`:
#'   no evaluation budget, so the run stops on the criterion or on
#'   `maxit`. Set a finite value to cap the cost of a run.
#' @param verbose Report progress? Defaults to `FALSE`.
#' @param refresh Report every this many iterations. Defaults to 100.
#' @param keep_trace Store the iteration path? Defaults to `FALSE`.
#'
#' @details
#' The idea is one line. Adam keeps an exponentially weighted average of the
#' gradient, \eqn{m_t}, and of its square, \eqn{v_t}, and steps
#' \deqn{x_{t+1} = x_t - \alpha\, \hat m_t / (\sqrt{\hat v_t} + \epsilon).}
#' The division is elementwise: a coordinate whose
#' gradient has been consistently large is divided by a large number and moves
#' modestly, while one whose gradient is small but persistent still moves. It is
#' a diagonal preconditioner assembled from the gradients already seen, so it
#' costs nothing beyond them. That is also its limit: a diagonal cannot
#' represent the correlation between parameters that [bfgs()] learns from the
#' same information.
#'
#' Both averages start at zero, so early on they are pulled towards it; dividing
#' by \eqn{1 - \beta^t} removes exactly that bias, and without it the first
#' iterations would barely move.
#'
#' # Relation to the descent methods
#'
#' Adam takes no line search and makes no attempt to decrease the objective
#' at every step. That freedom to go uphill is most of why it tolerates a
#' gradient that is right only on average, and it also means that none of the
#' usual reassurances apply: there is no guarantee of monotone progress, and
#' the run may end somewhere worse than it passed through.
#'
#' The practical consequence is that Adam is the wrong tool for a small
#' smooth problem where a Hessian is affordable. Measured on
#' \eqn{\sum(x - (1, 2))^2} from the origin, `adam(alpha = 0.1)` reaches the
#' answer in 2000 iterations and 2002 evaluations while [bfgs()] reaches it
#' in 2 and 3. Adam is for a long parameter vector, a noisy objective, or a
#' surface rough enough that a quadratic model is a fiction.
#'
#' # Stochastic objectives
#'
#' Adam draws no subsamples of its own, and that is the design. An optimizer
#' has no notion of an observation; a version that did would need a second
#' kind of objective to be told, a rule for which stopping rules such an
#' objective allows, and a way to report which was in force. A closure does
#' all of it in two lines, because a stochastic objective is an objective:
#'
#' ```r
#' batch <- function(par) {
#'   i <- sample.int(n, size = 0.05 * n)
#'   sum((y[i] - par)^2) / 2
#' }
#' minimize(adam(), batch, par = 0, gr = batch_gr)
#' ```
#'
#' Adam then behaves as it would on a minibatch of its own drawing, and
#' `set.seed()` governs the run because the draws happen in the caller's
#' code.
#'
#' **Resample inside the objective, not around the run.** Calling
#' `minimize(adam(maxit = 1), ...)` in a loop and drawing a new batch each
#' time does not work: \eqn{m} and \eqn{v} start at zero and the bias
#' correction restarts at \eqn{t = 1}, so every call takes a first step of
#' length \eqn{\alpha} and the accumulated moments, which are the whole of
#' the method, are thrown away at each one. Measured on the example below,
#' 200 one-iteration runs reach `2.90` where one 200-iteration run reaches
#' `2.956`, against a target of `2.986`.
#'
#' One thing to expect on such a run. The gradient consistency check
#' [minimize()] makes before starting compares `fn` and `gr` at one point,
#' and on a resampling objective those are two different minibatches, so the
#' check warns on every stochastic fit and the warning means nothing. Set
#' `options(optimizers7.check_gradient = FALSE)` around it.
#'
#' # Stopping
#'
#' The default criterion is [crit_never()], so the run ends when `maxit` is
#' reached and reports `converged = FALSE`, which is the truth about a run
#' nothing checked. With a fixed `alpha` Adam circles an optimum instead of
#' settling on it, so a tolerance on the gradient is usually a rule that
#' never fires.
#'
#' On an exact objective a real rule can be passed and will work. On a noisy
#' one nothing read from the objective or the gradient means much, both being
#' estimates, and the package cannot detect which kind it was given.
#'
#' # The safeguards
#'
#' `eps` floors the denominator, so a coordinate whose gradient has been
#' uniformly zero is not divided by zero.
#'
#' `decay` makes the learning rate \eqn{O(1/t)}, which is the Robbins-Monro
#' condition a run on a noisy objective needs to settle at the optimum
#' instead of rattling about it at a radius set by \eqn{\alpha}. In the
#' minibatch example below, five runs at `alpha = 0.05` scatter with a
#' standard deviation of `0.033` at `decay = 0` and `0.0079` at `0.01`. It
#' is off by default, an exact objective having nothing to average away.
#'
#' `amsgrad` replaces \eqn{v_t} by its running maximum. Reddi, Kale and Kumar
#' (2018) exhibited a convex problem on which Adam as published fails to
#' converge, because \eqn{v_t} can shrink and let a single large gradient
#' dominate the iterate long after it has passed. The maximum forbids that,
#' and the price is that a coordinate's effective rate
#' \eqn{\alpha/\sqrt{\max_s v_s}} can then never grow again. It is `FALSE` by
#' default so that `adam()` is Adam, and worth turning on when a run fails to
#' settle.
#'
#' A non-finite gradient or update ends the run rather than propagating a
#' `NaN` into every iterate after it. The result reports
#' `converged = FALSE` with the message `gradient not finite; stopped`, and
#' `par` is the last usable point.
#'
#' @return An S7 object of class [Adam], inheriting from [optimizer()], to be
#'   handed to [minimize()].
#'
#' @references
#' Kingma, D. P. and Ba, J. (2015). Adam: A Method for Stochastic Optimization.
#' *ICLR*.
#'
#' Reddi, S. J., Kale, S. and Kumar, S. (2018). On the Convergence of Adam and
#' Beyond. *ICLR*.
#'
#' @examples
#' adam()
#' adam(alpha = 0.05, amsgrad = TRUE)
#'
#' # On an exact quadratic it arrives, and pays for it: 2002 evaluations
#' # against BFGS's 3. The flag is FALSE because crit_never() checked nothing.
#' q <- function(p) sum((p - c(1, 2))^2)
#' qg <- function(p) 2 * (p - c(1, 2))
#' a <- minimize(adam(alpha = 0.1, maxit = 2000), q, c(0, 0), gr = qg)
#' c(a@par, evaluations = a@counts[["f"]], converged = a@converged)
#' minimize(bfgs(), q, c(0, 0), gr = qg)@counts[["f"]]
#'
#' # The objective it is for: the minibatch is drawn inside the function, so
#' # the optimizer never has to know about it. The check that compares fn
#' # against gr would compare two different batches, so it is turned off.
#' set.seed(1)
#' y <- rnorm(2000, mean = 3)
#' m <- 100
#' batch    <- function(p) { i <- sample.int(2000, m); sum((y[i] - p)^2) / 2 }
#' batch_gr <- function(p) { i <- sample.int(2000, m); -sum(y[i] - p) }
#'
#' old <- options(optimizers7.check_gradient = FALSE)
#' minimize(adam(alpha = 0.05, decay = 0.01, maxit = 2000),
#'          batch, par = 0, gr = batch_gr)@par
#' mean(y)
#'
#' # A decaying rate settles the run: the scatter over five runs falls with it.
#' spread <- function(d) sd(replicate(5,
#'   minimize(adam(alpha = 0.05, decay = d, maxit = 2000),
#'            batch, par = 0, gr = batch_gr)@par))
#' set.seed(2); c(none = spread(0), some = spread(0.01))
#' options(old)
#'
#' # A gradient that is not finite ends the run and says so.
#' bad <- minimize(adam(maxit = 50), q, c(0, 0), gr = function(p) c(NaN, 1))
#' c(bad@converged, bad@message)
#'
#' @seealso [bfgs()] for the smooth case, [crit_never()] for the default
#'   stopping rule, [sa()] for the other method that goes uphill on purpose.
#' @export
adam <- function(criterion = crit_never(),
                 alpha = 0.01, beta1 = 0.9, beta2 = 0.999, eps = 1e-8,
                 decay = 0, amsgrad = FALSE,
                 maxit = 1000, max_eval = Inf,
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

  Adam(
    name = "adam", criterion = criterion,
    maxit = maxit, max_eval = max_eval, verbose = verbose,
    refresh = refresh, keep_trace = keep_trace,
    alpha = alpha, beta1 = beta1, beta2 = beta2, eps = eps,
    decay = decay, amsgrad = amsgrad
  )
}


#' @title Minimize by Adam
#' @name minimize.Adam
#'
#' @description
#' Runs [adam()] on the objective: one gradient per iteration, the two
#' exponentially weighted moments updated from it, and a coordinatewise step
#' taken with no line search and no test that the objective fell.
#'
#' @param optimizer An `Adam` object.
#' @param fn,par,gr,he,lower,upper,... As in [minimize()]. `he` is accepted
#'   and ignored. Bounds are taken and removed by reparametrization.
#'
#' @return An [optimizer_result()]. Under the default [crit_never()] it
#'   reports `converged = FALSE` and `criterion_met` `iteration budget
#'   reached`, the run having ended on its budget.
#'
#' @keywords internal
S7::method(minimize, Adam) <-
  function(optimizer, fn, par, gr = NULL, he = NULL,
           lower = -Inf, upper = Inf, ...) {
    spec <- prepare_objective(optimizer, fn, par, gr, he)
    bounds <- check_bounds(lower, upper, par)

    t0 <- proc.time()[["elapsed"]]
    out <- adam_run(
      spec = spec, par = as.numeric(par),
      criterion = optimizer@criterion, crit_fn = crit_met,
      alpha = optimizer@alpha, beta1 = optimizer@beta1,
      beta2 = optimizer@beta2, eps = optimizer@eps,
      decay = optimizer@decay, amsgrad = optimizer@amsgrad,
      maxit = as.integer(optimizer@maxit),
      max_eval = budget_int(optimizer@max_eval),
      verbose = optimizer@verbose,
      refresh = as.integer(optimizer@refresh),
      keep_trace = optimizer@keep_trace,
      bounds = bounds
    )
    elapsed <- proc.time()[["elapsed"]] - t0

    build_result(out, optimizer, spec, elapsed)
  }
